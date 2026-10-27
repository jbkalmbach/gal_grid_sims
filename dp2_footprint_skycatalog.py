#!/usr/bin/env python
"""Cover the DP2 footprint with a hexagonal grid of exponential galaxies.

This builds a skyCatalog (for imSim) whose galaxies hex-tile the sky *inside*
each DP2 tract, so the union covers the real survey footprint without wasting
objects on the gaps between fields.  It reuses the per-object property handling
(scalar / scipy distribution / array) and flux computation from
``run_grid_skycatalog.py``.

The catalog is written by **streaming one tract at a time**: each tract is
tiled, its objects are split by healpix pixel, and each pixel's rows are
appended as a parquet row group.  Peak memory stays bounded to a single tract
(tens of thousands of objects) rather than the whole footprint (which can be
hundreds of millions of objects and was previously built entirely in memory).

Two entry points:

  * ``--save-tracts FILE`` queries the butler for every tract with DP2 data (or,
    with ``--field NAME``, only the tracts covered by that observing field) and
    writes a tract-list file, then exits.
  * ``--tracts-file FILE`` reads a tract-list file (full list or any subset) and
    builds the gridded skyCatalog covering those tracts.

Tract-list files are plain text: one ``tract`` per line (extra whitespace-
separated columns such as ra/dec centers are ignored), ``#`` comments allowed.

Environment (skyCatalogs, imSim, and the LSST stack for the butler):

    source ~/setup_aos_dev.sh
    setup -k -r /sdf/home/b/brycek/u/dev-repos/skyCatalogs
    setup -k -r /sdf/home/b/brycek/u/dev-repos/imSim

Example
-------
    # 1. Save the tract list for one field:
    python dp2_footprint_skycatalog.py --save-tracts field.txt \
        --field Rubin_SV_320_-15

    # 2. Build a 20-arcsec hex catalog covering those tracts:
    python dp2_footprint_skycatalog.py --tracts-file field.txt \
        --spacing 20 --mag "norm:loc=20,scale=0.5" --hlr 0.5 \
        --name dp2_sv320 --seed 42
"""
import argparse
import glob
import json
import os

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

# Reuse the building blocks that already back the notebook / grid script.
from run_grid_skycatalog import (
    N_TOPHAT, BANDS, resolve_param, mag0_band_fluxes,
    _galaxy_arrow_table, write_skycatalog_config,
)
from stars import (DEFAULT_STAR_DENSITY, DEFAULT_STAR_EXCLUSION_ARCSEC,
                   DEFAULT_STAR_MAG_MIN, DEFAULT_STAR_MAG_MAX, mag_spec)

# Defaults for the DP2 repo and the coadd collection we read the footprint from.
DP2_REPO = "dp2_prep"
DP2_COLLECTION = "LSSTCam/runs/DRP/DP2/v30_0_6_rc1/DM-53881/stage3"
DP2_SKYMAP = "lsst_cells_v2"
DP2_INSTRUMENT = "LSSTCam"


# --------------------------------------------------------------------------
# Butler footprint queries
# --------------------------------------------------------------------------
def _open_skymap(butler):
    return butler.get("skyMap", skymap=DP2_SKYMAP, collections="skymaps")


def query_footprint_tracts(field=None, repo=DP2_REPO, collection=DP2_COLLECTION):
    """Return a sorted list of DP2 tracts that have coadd (object_patch) data.

    If `field` is given (a visit ``target_name``, e.g. ``Rubin_SV_320_-15``),
    restrict to tracts overlapping that field's visits.
    """
    from lsst.daf.butler import Butler

    butler = Butler(repo)
    reg = butler.registry

    refs = reg.queryDatasets("object_patch", collections=collection,
                             findFirst=True)
    tracts_with_data = {r.dataId["tract"] for r in refs}

    if field is None:
        return sorted(tracts_with_data)

    visits = [v.id for v in reg.queryDimensionRecords(
        "visit", instrument=DP2_INSTRUMENT,
        where="visit.target_name = :nm", bind={"nm": field})]
    if not visits:
        raise ValueError(f"no visits found with target_name = {field!r}")
    field_tracts = set()
    for v in visits:
        for did in reg.queryDataIds(["tract"], skymap=DP2_SKYMAP,
                                    instrument=DP2_INSTRUMENT, visit=v):
            field_tracts.add(did["tract"])
    return sorted(field_tracts & tracts_with_data)


def save_tract_list(path, tracts, repo=DP2_REPO, collection=DP2_COLLECTION,
                    field=None):
    """Write a tract-list file with tract id and sky center for each tract."""
    from lsst.daf.butler import Butler

    butler = Butler(repo)
    skymap = _open_skymap(butler)
    with open(path, "w") as f:
        f.write("# DP2 tracts with object_patch data\n")
        if field is not None:
            f.write(f"# field (visit.target_name): {field}\n")
        f.write(f"# collection: {collection}\n")
        f.write(f"# skymap: {DP2_SKYMAP}\n")
        f.write("# tract  ra_center_deg  dec_center_deg\n")
        for t in tracts:
            c = skymap[t].getCtrCoord()
            f.write(f"{t} {c.getRa().asDegrees():.5f} "
                    f"{c.getDec().asDegrees():.5f}\n")
    return path


def read_tract_list(path):
    """Read tract ids from a tract-list file (first column, '#' comments)."""
    tracts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tracts.append(int(line.split()[0]))
    return sorted(set(tracts))


# --------------------------------------------------------------------------
# Hex tiling of a tract's inner sky region
# --------------------------------------------------------------------------
def hex_tile_tract(tract_info, spacing_arcsec):
    """Return (ra, dec) arrays [deg] hex-tiling a tract's inner sky region.

    Uses the tract's non-overlapping inner region so adjacent tracts tile the
    sky without double-counting.  Alternate rows are offset by half a step and
    rows are spaced by sqrt(3)/2 * step, giving nearest-neighbor separation
    `spacing_arcsec` (in the local tangent plane).
    """
    import lsst.sphgeom as sph

    region = tract_info.inner_sky_region  # sphgeom Box
    bb = region.getBoundingBox()
    ra_c = bb.getCenter().getLon().asDegrees()
    dec_c = bb.getCenter().getLat().asDegrees()
    ra_w = bb.getWidth().asDegrees()
    dec_h = bb.getHeight().asDegrees()

    step = spacing_arcsec / 3600.0  # deg
    dec_lo, dec_hi = dec_c - dec_h / 2.0, dec_c + dec_h / 2.0
    row_decs = np.arange(dec_lo, dec_hi + step, step * np.sqrt(3.0) / 2.0)

    ras, decs = [], []
    for j, dd in enumerate(row_decs):
        cosd = np.cos(np.deg2rad(dd))
        rstep = step / cosd                    # keep on-sky spacing uniform
        off = 0.5 * rstep * (j % 2)            # offset odd rows -> hexagonal
        ra_lo = ra_c - ra_w / 2.0 / cosd
        ra_hi = ra_c + ra_w / 2.0 / cosd
        row_ras = np.arange(ra_lo + off, ra_hi + rstep, rstep)
        for r in row_ras:
            uv = sph.UnitVector3d(sph.LonLat.fromDegrees(r % 360.0, dd))
            if region.contains(uv):
                ras.append(r % 360.0)
                decs.append(dd)
    return np.array(ras), np.array(decs)


# galaxy_id = tract * ID_STRIDE + local index.  Order-independent and stable
# per tract, so ids stay globally unique across resumes without bookkeeping.
# (A tract's object count must stay below this; ~10M objects/tract needs a
# spacing near 2 arcsec, far finer than any practical grid.)
ID_STRIDE = 100_000_000


# --------------------------------------------------------------------------
# Build the per-object chunk for a single tract
# --------------------------------------------------------------------------
def _tract_chunk(ra, dec, rng, id_base, *, mag, hlr, ellipticity, pa_deg,
                 g1, g2, redshift):
    """Build the per-object column dict for one tract's (ra, dec) points.

    Every non-position property is resolved through ``resolve_param`` (scalar /
    scipy distribution / array) using this tract's `rng`.  `galaxy_id` runs from
    `id_base` so it is globally unique and stable regardless of tract order.
    """
    n = ra.size
    if n >= ID_STRIDE:
        raise ValueError(f"tract has {n} objects, exceeds ID_STRIDE={ID_STRIDE}; "
                         "use coarser spacing or raise ID_STRIDE")
    mag_v = resolve_param(mag, n, rng)
    hlr_v = resolve_param(hlr, n, rng)
    e_v = resolve_param(ellipticity, n, rng)
    pa_v = np.deg2rad(resolve_param(pa_deg, n, rng))
    g1_v = resolve_param(g1, n, rng)
    g2_v = resolve_param(g2, n, rng)
    z_v = resolve_param(redshift, n, rng)

    e1 = e_v * np.cos(2.0 * pa_v)
    e2 = e_v * np.sin(2.0 * pa_v)

    return {
        "galaxy_id": np.arange(id_base, id_base + n, dtype=np.int64),
        "ra": ra.astype(np.float64),
        "dec": dec.astype(np.float64),
        "redshift": z_v.astype(np.float64),
        "redshift_hubble": z_v.astype(np.float64),
        "shear_1": g1_v.astype(np.float64),
        "shear_2": g2_v.astype(np.float64),
        "convergence": np.zeros(n, dtype=np.float64),
        "size_disk_true": hlr_v.astype(np.float32),
        "size_minor_disk_true": hlr_v.astype(np.float32),
        "sersic_disk": np.ones(n, dtype=np.float32),
        "ellipticity_1_disk_true": e1.astype(np.float64),
        "ellipticity_2_disk_true": e2.astype(np.float64),
        # sed_val_disk is a constant flat vector; build lazily per row group.
        "MW_av": np.zeros(n, dtype=np.float32),
        "MW_rv": np.full(n, 3.1, dtype=np.float32),
        "mag": mag_v.astype(np.float64),
    }


def _chunk_arrow_tables(chunk, idx, f0, bands):
    """Return (galaxy_table, flux_table) pyarrow tables for chunk rows `idx`."""
    # sed_val_disk (flat tophat, shape only) built on the fly for these rows.
    chunk = dict(chunk)
    chunk["sed_val_disk"] = [list(np.ones(N_TOPHAT))] * len(chunk["ra"])
    gal = _galaxy_arrow_table(chunk, idx)

    f32 = pa.float32()
    mag = chunk["mag"][idx]
    cols = {"galaxy_id": pa.array(chunk["galaxy_id"][idx], type=pa.int64())}
    for b in bands:
        fb = (f0[b] * 10.0 ** (-0.4 * mag)).astype(np.float32)
        cols[f"lsst_flux_{b}"] = pa.array(fb, type=f32)
    return gal, pa.table(cols)


# --------------------------------------------------------------------------
# Resume manifest
# --------------------------------------------------------------------------
def _manifest_path(cat_dir):
    return os.path.join(cat_dir, "_progress.json")


def _load_manifest(cat_dir):
    """Return the resume manifest dict, or a fresh one if none / incompatible."""
    path = _manifest_path(cat_dir)
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {"done_tracts": [], "n_total": 0, "params": None}


def _save_manifest(cat_dir, manifest):
    """Atomically write the manifest (temp file + rename)."""
    path = _manifest_path(cat_dir)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(manifest, f)
    os.replace(tmp, path)


def _run_params(spacing_arcsec, nside, bands, seed, star_density, star_seed,
                star_exclusion_arcsec):
    """Signature of run parameters; a resume must match these."""
    return {"spacing_arcsec": spacing_arcsec, "nside": nside,
            "bands": list(bands), "seed": seed,
            "star_density": star_density, "star_seed": star_seed,
            "star_exclusion_arcsec": star_exclusion_arcsec}


# --------------------------------------------------------------------------
# Stream the footprint to per-tract shards, then merge into per-healpix files
# --------------------------------------------------------------------------
def _write_tract_shards(tract, ra, dec, chunk, f0, bands, nside, shard_dir):
    """Write one tract's objects as per-healpix shard files; return pixel set.

    Each shard (``<hp>/galaxy_tract<t>.parquet`` and the flux counterpart) is a
    self-contained parquet file closed on write, so a kill leaves only complete
    shards -- nothing to corrupt or truncate.
    """
    import healpy as hp
    pix = hp.ang2pix(nside, ra, dec, lonlat=True, nest=False)
    pixels = np.unique(pix)
    for hpix in pixels:
        idx = np.where(pix == hpix)[0]
        gal, flux = _chunk_arrow_tables(chunk, idx, f0, bands)
        pdir = os.path.join(shard_dir, str(int(hpix)))
        os.makedirs(pdir, exist_ok=True)
        pq.write_table(gal, os.path.join(pdir, f"galaxy_tract{tract}.parquet"))
        pq.write_table(flux, os.path.join(pdir, f"galaxy_flux_tract{tract}.parquet"))
    return [int(p) for p in pixels]


def _write_star_shards(tract, ids, ra, dec, mag, f0, bands, nside, shard_dir):
    """Write one tract's stars as per-healpix pointsource shard files.

    Mirrors ``_write_tract_shards`` for the point-source population: a
    ``<hp>/pointsource_tract<t>.parquet`` main file and its flux counterpart,
    built from the SAME per-pixel index so skyCatalogs' positional main/flux
    join stays aligned.  Returns the pixel set written (empty if no stars).
    """
    import healpy as hp
    from stars import star_arrow_tables

    if ids.size == 0:
        return []
    pix = hp.ang2pix(nside, ra, dec, lonlat=True, nest=False)
    pixels = np.unique(pix)
    for hpix in pixels:
        idx = np.where(pix == hpix)[0]
        main, flux = star_arrow_tables(ids[idx], ra[idx], dec[idx], mag[idx],
                                       f0, bands)
        pdir = os.path.join(shard_dir, str(int(hpix)))
        os.makedirs(pdir, exist_ok=True)
        pq.write_table(main, os.path.join(pdir, f"pointsource_tract{tract}.parquet"))
        pq.write_table(flux, os.path.join(pdir, f"pointsource_flux_tract{tract}.parquet"))
    return [int(p) for p in pixels]


def _merge_shards(cat_dir, shard_dir, delete=True, logger=print,
                  kinds=("galaxy", "galaxy_flux")):
    """Merge per-tract shards into one <kind>_<hp>.parquet per healpix pixel.

    Merges pixel by pixel and, when `delete` is True, removes each shard right
    after it is read and the pixel's shard directory once merged.  This bounds
    the transient extra disk to a single pixel's data rather than a full second
    copy of the whole footprint.

    `kinds` lists the shard prefixes to merge (galaxy + its flux, and, when
    stars are present, pointsource + its flux).  For each kind the shards are
    merged in the same sorted (tract) order so the star main/flux files stay
    positionally aligned for skyCatalogs' join.
    """
    if not os.path.isdir(shard_dir):
        return
    pixels = sorted(int(d) for d in os.listdir(shard_dir)
                    if d.isdigit() and os.path.isdir(os.path.join(shard_dir, d)))
    if logger:
        logger(f"merging shards for {len(pixels)} healpix pixels"
               + (" (deleting as we go)" if delete else "") + "...")
    for hpix in pixels:
        pdir = os.path.join(shard_dir, str(hpix))
        for kind in kinds:
            shards = sorted(glob.glob(os.path.join(pdir, f"{kind}_tract*.parquet")))
            if not shards:
                continue
            out = os.path.join(cat_dir, f"{kind}_{hpix}.parquet")
            writer = None
            try:
                for s in shards:
                    tbl = pq.read_table(s)
                    if writer is None:
                        writer = pq.ParquetWriter(out, tbl.schema)
                    writer.write_table(tbl)
                    if delete:
                        os.remove(s)   # free each shard as soon as it's merged
            finally:
                if writer is not None:
                    writer.close()
        if delete:
            # Pixel fully merged; drop its (now-empty) shard directory.
            try:
                os.rmdir(pdir)
            except OSError:
                pass


def write_footprint_skycatalog(tracts, spacing_arcsec, out_dir, *,
                               catalog_name="dp2_footprint", repo=DP2_REPO,
                               nside=32, bands=BANDS, mag=20.0, hlr=0.5,
                               ellipticity=0.0, pa_deg=0.0, g1=0.0, g2=0.0,
                               redshift=0.1, seed=None, star_density=0.0,
                               star_mag=22.0, star_exclusion_arcsec=0.0,
                               star_seed=None, resume=False,
                               keep_shards=False, logger=print):
    """Hex-tile each tract and stream the objects to parquet, resumably.

    Each tract is tiled, resolved, split by healpix pixel, and written as
    self-contained per-tract shard files (so a killed job leaves only complete
    shards).  Completed tracts are recorded in a manifest; with ``resume=True``
    those are skipped on a rerun.  When all tracts are done, shards are merged
    into one ``galaxy_<hp>.parquet`` / ``galaxy_flux_<hp>.parquet`` per pixel and
    the skyCatalog.yaml is written.  Returns (config_path, n_objects_total).

    If `star_density` > 0 (stars per square arcminute), each tract also gets a
    sparse, irregular (uniform-random) population of point sources written as
    ``pointsource_<hp>.parquet`` / ``pointsource_flux_<hp>.parquet``, and a
    ``star`` block is added to the config.  These stars are what let the LSST
    single-frame calibration (PSF modeling + ``matchPessimisticB`` astrometry)
    run on an otherwise degenerate galaxy lattice.  `star_mag` is a scalar or
    scipy spec; `star_seed` (default: `seed`) drives a star RNG stream that is
    independent of the galaxy draws.  If `star_exclusion_arcsec` > 0, candidate
    stars within that angular radius of a grid galaxy in the tract are rejected
    so stars stay isolated point sources (keep it well below half `spacing`).

    Memory stays bounded to a single tract.  Per-tract RNG seeding
    (``seed`` combined with the tract id) and ``galaxy_id = tract*ID_STRIDE + i``
    make the result independent of tract order, so a resumed run reproduces the
    same catalog as an uninterrupted one.
    """
    import time
    from lsst.daf.butler import Butler
    from stars import (n_stars_for_area, draw_positions_in_region,
                       STAR_ID_NAMESPACE, STAR_SEED_TAG)

    with_stars = star_density and star_density > 0
    if star_seed is None:
        star_seed = seed

    cat_dir = os.path.join(out_dir, catalog_name)
    shard_dir = os.path.join(cat_dir, "_shards")
    os.makedirs(shard_dir, exist_ok=True)

    params = _run_params(spacing_arcsec, nside, bands, seed,
                         star_density if with_stars else 0.0,
                         star_seed if with_stars else None,
                         star_exclusion_arcsec if with_stars else 0.0)
    manifest = _load_manifest(cat_dir) if resume else \
        {"done_tracts": [], "n_total": 0, "n_stars_total": 0, "params": params}
    if resume and manifest["params"] not in (None, params):
        raise ValueError(
            "resume parameters do not match the existing run:\n"
            f"  existing: {manifest['params']}\n  requested: {params}\n"
            "Use a fresh --name or matching parameters.")
    manifest["params"] = params
    manifest.setdefault("n_stars_total", 0)
    done = set(manifest["done_tracts"])
    if resume and done and logger:
        logger(f"resuming: {len(done)} of {len(tracts)} tracts already done")

    butler = Butler(repo)
    skymap = _open_skymap(butler)
    f0 = mag0_band_fluxes(bands)          # compute SED normalizations once

    todo = [t for t in tracts if t not in done]
    n_total = manifest["n_total"]
    n_stars_total = manifest["n_stars_total"]
    t_start = time.time()
    for i, t in enumerate(todo, 1):
        ra, dec = hex_tile_tract(skymap[t], spacing_arcsec)
        if ra.size:
            # Per-tract seed so results are order-independent and resumable.
            tract_rng = np.random.default_rng(None if seed is None
                                              else (seed, t))
            id_base = t * ID_STRIDE
            chunk = _tract_chunk(ra, dec, tract_rng, id_base, mag=mag, hlr=hlr,
                                 ellipticity=ellipticity, pa_deg=pa_deg,
                                 g1=g1, g2=g2, redshift=redshift)
            npix = len(_write_tract_shards(t, ra, dec, chunk, f0, bands,
                                           nside, shard_dir))
            n_total += ra.size
        else:
            npix = 0

        # Draw this tract's sparse irregular stars from a separate RNG stream so
        # the star realization does not depend on the galaxy param draw order.
        n_stars = 0
        if with_stars:
            region = skymap[t].inner_sky_region
            area_sq_deg = region.getArea() * (180.0 / np.pi) ** 2  # sr -> deg^2
            n_want = n_stars_for_area(area_sq_deg, star_density)
            if n_want > 0:
                star_rng = np.random.default_rng(
                    None if star_seed is None
                    else (star_seed, t, STAR_SEED_TAG))
                # Exclude stars landing on this tract's grid galaxies (ra, dec).
                sra, sdec = draw_positions_in_region(
                    region, n_want, star_rng, galaxy_ra=ra, galaxy_dec=dec,
                    exclusion_arcsec=star_exclusion_arcsec)
                if sra.size:
                    smag = resolve_param(star_mag, sra.size, star_rng)
                    sids = STAR_ID_NAMESPACE + t * ID_STRIDE + \
                        np.arange(sra.size, dtype=np.int64)
                    _write_star_shards(t, sids, sra, sdec, smag, f0, bands,
                                       nside, shard_dir)
                    n_stars = sra.size
                    n_stars_total += n_stars

        # Record progress after each tract so a kill loses at most one tract.
        manifest["done_tracts"].append(t)
        manifest["n_total"] = n_total
        manifest["n_stars_total"] = n_stars_total
        _save_manifest(cat_dir, manifest)

        if logger:
            elapsed = time.time() - t_start
            rate = i / elapsed if elapsed > 0 else 0.0
            eta = (len(todo) - i) / rate if rate > 0 else float("nan")
            star_msg = f", {n_stars} stars" if with_stars else ""
            logger(f"[{i}/{len(todo)}] tract {t}: {ra.size} objects{star_msg} "
                   f"-> {npix} pixels | total {n_total} gal / "
                   f"{n_stars_total} stars | "
                   f"{rate*60:.1f} tracts/min | ETA {eta/60:.1f} min")

    if logger:
        logger(f"all {len(tracts)} tracts done ({n_total} galaxies, "
               f"{n_stars_total} stars); finalizing")
    # Merge pixel by pixel, deleting each shard as it is consumed (unless the
    # user wants to keep them), so the transient extra disk never exceeds one
    # pixel's worth rather than a full second copy of the footprint.
    kinds = ["galaxy", "galaxy_flux"]
    if with_stars:
        kinds += ["pointsource", "pointsource_flux"]
    _merge_shards(cat_dir, shard_dir, delete=not keep_shards, logger=logger,
                  kinds=kinds)
    config_path = write_skycatalog_config(cat_dir, nside=nside,
                                          with_stars=bool(with_stars))

    if not keep_shards:
        # _merge_shards already removed the per-pixel shards + dirs; clear the
        # now-empty shard root and the manifest.
        import shutil
        shutil.rmtree(shard_dir, ignore_errors=True)
        mpath = _manifest_path(cat_dir)
        if os.path.exists(mpath):
            os.remove(mpath)
        if logger:
            logger("removed shard files and manifest (use --keep-shards to keep)")
    return config_path, n_total, n_stars_total


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def _param(s):
    """Parse a CLI value that is a float or a '<dist>:k=v,...' scipy spec."""
    from run_grid_skycatalog import parse_param_spec
    return parse_param_spec(s)


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    # Mode 1: save a tract list and exit.
    p.add_argument("--save-tracts", metavar="FILE",
                   help="query the butler, write a tract-list file, and exit")
    p.add_argument("--field", default=None,
                   help="restrict --save-tracts to a visit target_name field, "
                        "e.g. Rubin_SV_320_-15 (default: all DP2 tracts)")
    # Mode 2: build a catalog from a tract list.
    p.add_argument("--tracts-file", metavar="FILE",
                   help="tract-list file to build the gridded catalog from")
    p.add_argument("--repo", default=DP2_REPO)
    p.add_argument("--collection", default=DP2_COLLECTION)
    # Grid + galaxy properties.
    p.add_argument("--spacing", type=float, default=20.0,
                   help="hex nearest-neighbor spacing [arcsec] (default 20)")
    p.add_argument("--mag", type=_param, default=20.0)
    p.add_argument("--hlr", type=_param, default=0.5)
    p.add_argument("--ellipticity", type=_param, default=0.0)
    p.add_argument("--pa-deg", type=_param, default=0.0)
    p.add_argument("--g1", type=_param, default=0.0)
    p.add_argument("--g2", type=_param, default=0.0)
    p.add_argument("--redshift", type=_param, default=0.1)
    p.add_argument("--seed", type=int, default=None)
    # Sparse irregular star population (for LSST calibration).  On by default at
    # the module default density; pass --star-density 0 for galaxies only.
    p.add_argument("--star-density", type=float, default=None,
                   help="point-source surface density [stars/arcmin^2]; adds a "
                        "sparse irregular star population for calibration "
                        f"(default: {DEFAULT_STAR_DENSITY}; 0 = galaxies only)")
    p.add_argument("--star-mag-min", type=float, default=DEFAULT_STAR_MAG_MIN,
                   help="bright end of the uniform star AB magnitude range "
                        f"(default {DEFAULT_STAR_MAG_MIN})")
    p.add_argument("--star-mag-max", type=float, default=DEFAULT_STAR_MAG_MAX,
                   help="faint end of the uniform star AB magnitude range "
                        f"(default {DEFAULT_STAR_MAG_MAX}); set equal to "
                        "--star-mag-min for a single fixed magnitude")
    p.add_argument("--star-mag", type=_param, default=None,
                   help="override the min/max range with an explicit star AB "
                        "magnitude: float or scipy dist spec (advanced)")
    p.add_argument("--star-exclusion-arcsec", type=float,
                   default=DEFAULT_STAR_EXCLUSION_ARCSEC,
                   help="reject stars within this angular radius of a grid "
                        f"galaxy (default {DEFAULT_STAR_EXCLUSION_ARCSEC}; 0 "
                        "disables). Keep well below half the --spacing.")
    p.add_argument("--star-seed", type=int, default=None,
                   help="RNG seed for star positions/mags (default: --seed)")
    p.add_argument("--nside", type=int, default=32)
    p.add_argument("--out-dir",
                   default="/sdf/data/rubin/user/brycek/grid_galaxies/catalogs")
    p.add_argument("--name", default="dp2_footprint")
    p.add_argument("--resume", action="store_true",
                   help="skip tracts already recorded in the run's manifest "
                        "(_progress.json) and continue from where a prior run "
                        "left off")
    p.add_argument("--keep-shards", action="store_true",
                   help="keep the per-tract shard files and manifest after the "
                        "final merge (default: remove them to save disk)")
    args = p.parse_args()

    if args.save_tracts:
        tracts = query_footprint_tracts(field=args.field, repo=args.repo,
                                        collection=args.collection)
        save_tract_list(args.save_tracts, tracts, repo=args.repo,
                        collection=args.collection, field=args.field)
        label = f"field {args.field}" if args.field else "full DP2 footprint"
        print(f"Wrote {len(tracts)} tracts ({label}) to {args.save_tracts}")
        return

    if not args.tracts_file:
        p.error("give either --save-tracts (to write a list) or "
                "--tracts-file (to build a catalog)")

    star_density = (DEFAULT_STAR_DENSITY if args.star_density is None
                    else args.star_density)
    # Explicit --star-mag (if given) wins; otherwise draw uniformly across
    # the [--star-mag-min, --star-mag-max] range.
    star_mag = (args.star_mag if args.star_mag is not None
                else mag_spec(args.star_mag_min, args.star_mag_max))
    tracts = read_tract_list(args.tracts_file)
    print(f"Building grid over {len(tracts)} tracts "
          f"at {args.spacing} arcsec spacing (streaming, one tract at a time)...")
    if star_density > 0:
        mag_desc = (f"mag {args.star_mag}" if args.star_mag is not None
                    else f"mag {args.star_mag_min}-{args.star_mag_max}")
        print(f"Adding sparse irregular stars at {star_density} stars/arcmin^2 "
              f"({mag_desc}), excluding within "
              f"{args.star_exclusion_arcsec}\" of a galaxy.")
    config_path, n_total, n_stars_total = write_footprint_skycatalog(
        tracts, args.spacing, args.out_dir, catalog_name=args.name,
        repo=args.repo, nside=args.nside,
        mag=args.mag, hlr=args.hlr, ellipticity=args.ellipticity,
        pa_deg=args.pa_deg, g1=args.g1, g2=args.g2, redshift=args.redshift,
        seed=args.seed, star_density=star_density, star_mag=star_mag,
        star_exclusion_arcsec=args.star_exclusion_arcsec,
        star_seed=args.star_seed, resume=args.resume,
        keep_shards=args.keep_shards)
    print(f"Wrote {n_total} galaxies and {n_stars_total} stars "
          f"over {len(tracts)} tracts.")
    print(f"skyCatalog config: {config_path}")


if __name__ == "__main__":
    main()

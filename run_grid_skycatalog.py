#!/usr/bin/env python
"""Build a grid of simple exponential galaxies as a skyCatalog for imSim.

Each galaxy is a single exponential (Sersic n=1) disk. Per-object properties
(magnitude, half-light radius, ellipticity, position angle, redshift, and the
applied shear) may each be supplied as:

  * a single scalar          -> the same value for every galaxy, or
  * a scipy.stats frozen distribution (anything with an ``.rvs`` method)
                             -> an independent random draw per galaxy, or
  * an explicit array of length ``n*n`` -> per-object values in object order.

This lets you, e.g., fix the half-light radius but draw magnitudes from
``scipy.stats.norm(loc=20, scale=0.5)``.  All random draws share one seeded
``numpy`` Generator so a given ``seed`` reproduces the catalog exactly.

A constant (or drawn) applied shear (g1, g2) is injected via the skyCatalogs
``shear_1``/``shear_2`` columns for shear-recovery tests.  Fluxes are written
for all ugrizy bands.

The catalog is consumed unmodified by imSim's `sky_catalog` input
(imsim/skycat.py), which only relies on the skyCatalogs GalaxyObject methods
`get_gsobject_components`, `get_observer_sed_components`, and `get_LSST_flux`.

Environment (skyCatalogs & imSim are not in the default stack setup):

    source ~/setup_aos_dev.sh
    setup -k -r /sdf/home/b/brycek/u/dev-repos/skyCatalogs
    setup -k -r /sdf/home/b/brycek/u/dev-repos/imSim

This module mirrors, cell-for-cell, the code in ``grid_skycatalog_demo.ipynb``.
"""
import argparse
import os

import numpy as np
import scipy.stats
import yaml
import healpy as hp
import pyarrow as pa
import pyarrow.parquet as pq
import galsim

from skycatalogs.objects.base_object import load_lsst_bandpasses

# --------------------------------------------------------------------------
# Fixed catalog metadata
# --------------------------------------------------------------------------

# cosmoDC2-style 30-bin tophat SED definition (start, width) in angstroms.
# Taken from skyCatalogs ci_sample/skyCatalog.yaml.  The GalaxyObject SED path
# requires this; for our galaxies the SED is flat (all bins equal) and only
# sets the *shape* -- the total flux is taken from the lsst_flux_* columns.
TOPHAT_BINS = [
    [1000, 246], [1246, 306], [1552, 381], [1933, 474], [2407, 591],
    [2998, 186], [3184, 197], [3381, 209], [3590, 222], [3812, 236],
    [4048, 251], [4299, 266], [4565, 283], [4848, 300], [5148, 319],
    [5467, 339], [5806, 360], [6166, 382], [6548, 406], [6954, 431],
    [7385, 458], [7843, 486], [8329, 517], [8846, 549], [9395, 583],
    [9978, 1489], [11467, 1710], [13177, 1966], [15143, 2259], [17402, 2596],
]
N_TOPHAT = len(TOPHAT_BINS)

COSMOLOGY = {"H0": 71.0, "Ob0": 0.0448, "Om0": 0.2648, "n_s": 0.963, "sigma8": 0.8}

BANDS = "ugrizy"


# --------------------------------------------------------------------------
# Parameter resolution: scalar | scipy distribution | per-object array
# --------------------------------------------------------------------------
def resolve_param(spec, n, rng):
    """Turn a parameter `spec` into a length-`n` array of per-object values.

    Parameters
    ----------
    spec :
        One of: a scalar; a scipy.stats frozen distribution (any object with an
        ``.rvs`` method); or an array-like of length `n`.
    n : int
        Number of objects.
    rng : numpy.random.Generator
        Source of randomness for distribution draws (via ``random_state=rng``).

    Returns
    -------
    numpy.ndarray of float, shape (n,)
    """
    if hasattr(spec, "rvs"):  # scipy.stats frozen distribution
        return np.asarray(spec.rvs(size=n, random_state=rng), dtype=float)
    arr = np.asarray(spec, dtype=float)
    if arr.ndim == 0:
        return np.full(n, float(arr))
    if arr.shape != (n,):
        raise ValueError(
            f"array parameter has shape {arr.shape}, expected ({n},)")
    return arr


def parse_param_spec(s):
    """Parse a CLI parameter that is either a float or a scipy distribution.

    Distribution syntax: ``"<distname>:key=val,key=val"``, e.g.
    ``"norm:loc=20,scale=0.5"`` or ``"uniform:loc=0.2,scale=0.8"``.  `distname`
    is any distribution in ``scipy.stats``; keyword values are parsed as floats.
    """
    try:
        return float(s)
    except ValueError:
        pass
    if ":" not in s:
        raise argparse.ArgumentTypeError(
            f"{s!r} is neither a float nor a '<dist>:k=v,...' spec")
    name, _, kvs = s.partition(":")
    try:
        dist = getattr(scipy.stats, name.strip())
    except AttributeError:
        raise argparse.ArgumentTypeError(
            f"unknown scipy.stats distribution {name!r}")
    kwargs = {}
    for tok in kvs.split(","):
        if not tok.strip():
            continue
        k, _, v = tok.partition("=")
        kwargs[k.strip()] = float(v)
    return dist(**kwargs)


# --------------------------------------------------------------------------
# 1. Build the grid of galaxy parameters
# --------------------------------------------------------------------------
def build_grid(ra0, dec0, n, spacing_arcsec, *,
               mag=20.0, hlr=0.5, ellipticity=0.0, pa_deg=0.0,
               g1=0.0, g2=0.0, redshift=0.1, seed=None):
    """Return a dict of per-object column arrays for an n x n galaxy grid.

    Galaxies sit on a hexagonal n x n grid centered on (`ra0`, `dec0`) with
    nearest-neighbor separation `spacing_arcsec` (alternate rows offset by half
    a step, rows spaced by sqrt(3)/2 * step).  Each of `mag`, `hlr`, `ellipticity`,
    `pa_deg`, `g1`, `g2`, and `redshift` may be a scalar, a scipy.stats frozen
    distribution, or a length-`n*n` array (see `resolve_param`).

    Parameters
    ----------
    ra0, dec0 : float
        Grid center on the sky, in degrees.
    n : int
        Number of galaxies per side (total = n*n).
    spacing_arcsec : float
        On-sky spacing between adjacent grid points, in arcsec.
    mag : scalar | scipy dist | array
        AB magnitude (same in all bands -> flat SED color).
    hlr : scalar | scipy dist | array
        Half-light radius, arcsec.
    ellipticity : scalar | scipy dist | array
        Intrinsic ellipticity magnitude |e|.
    pa_deg : scalar | scipy dist | array
        Position angle of the intrinsic ellipticity, degrees.
    g1, g2 : scalar | scipy dist | array
        Applied (cosmological) shear injected into each galaxy.
    redshift : scalar | scipy dist | array
        Galaxy redshift.  Must be > 0 (the tophat SED divides by the luminosity
        distance).  The SED is flat and its normalization is overridden by the
        lsst_flux_* columns, so the value only mildly affects the chromatic
        rendering.
    seed : int or None
        Seed for the numpy Generator used for all distribution draws.
    """
    rng = np.random.default_rng(seed)
    n_obj = n * n

    # Hexagonal grid: n x n points where alternate rows are offset by half a
    # step in x and rows are spaced by sqrt(3)/2 * step in y.  With this layout
    # the nearest-neighbor distance equals `spacing_arcsec` exactly.
    step = spacing_arcsec / 3600.0  # deg
    ix, iy = np.meshgrid(np.arange(n), np.arange(n), indexing="xy")
    ix = ix.ravel()
    iy = iy.ravel()

    x = ix.astype(float) * step + 0.5 * step * (iy % 2)  # offset odd rows
    y = iy.astype(float) * step * (np.sqrt(3.0) / 2.0)
    x -= x.mean()  # center the grid on (ra0, dec0)
    y -= y.mean()

    # RA offset scaled by cos(dec) so on-sky spacing is uniform.
    ra = ra0 + x / np.cos(np.deg2rad(dec0))
    dec = dec0 + y

    # Resolve per-object parameters.  Order is fixed so a given seed is
    # reproducible regardless of which parameters are distributions.
    mag_v = resolve_param(mag, n_obj, rng)
    hlr_v = resolve_param(hlr, n_obj, rng)
    e_v = resolve_param(ellipticity, n_obj, rng)
    pa_v = np.deg2rad(resolve_param(pa_deg, n_obj, rng))
    g1_v = resolve_param(g1, n_obj, rng)
    g2_v = resolve_param(g2, n_obj, rng)
    z_v = resolve_param(redshift, n_obj, rng)

    e1 = e_v * np.cos(2.0 * pa_v)
    e2 = e_v * np.sin(2.0 * pa_v)

    # Flat tophat SED (shape only; flux comes from lsst_flux_* columns).
    sed_val = [list(np.ones(N_TOPHAT))] * n_obj

    table = {
        "galaxy_id": np.arange(n_obj, dtype=np.int64),
        "ra": ra.astype(np.float64),
        "dec": dec.astype(np.float64),
        "redshift": z_v.astype(np.float64),
        "redshift_hubble": z_v.astype(np.float64),
        "shear_1": g1_v.astype(np.float64),
        "shear_2": g2_v.astype(np.float64),
        "convergence": np.zeros(n_obj, dtype=np.float64),
        "size_disk_true": hlr_v.astype(np.float32),
        "size_minor_disk_true": hlr_v.astype(np.float32),
        "sersic_disk": np.ones(n_obj, dtype=np.float32),
        "ellipticity_1_disk_true": e1.astype(np.float64),
        "ellipticity_2_disk_true": e2.astype(np.float64),
        "sed_val_disk": sed_val,
        "MW_av": np.zeros(n_obj, dtype=np.float32),
        "MW_rv": np.full(n_obj, 3.1, dtype=np.float32),
        # Book-keeping (not read by skyCatalogs, handy for analysis):
        "mag": mag_v.astype(np.float64),
    }
    return table


# --------------------------------------------------------------------------
# 2. Magnitude -> LSST photon fluxes
# --------------------------------------------------------------------------
def compute_lsst_fluxes(mag, bands=BANDS):
    """Return {band: flux_array} in photons/s/cm^2 for AB source(s) of `mag`.

    `mag` may be a scalar or an array; the return values are arrays with the
    same length (a scalar is treated as length 1).  Because the SED shape is
    fixed, flux scales exactly as ``10**(-0.4*mag)`` from the mag=0 flux, so
    only one SED integration per band is needed regardless of the number of
    objects.  Uses the same LSST bandpasses skyCatalogs uses, so the values are
    on the exact footing that `get_LSST_flux` / imSim expect.
    """
    f0 = mag0_band_fluxes(bands)
    mag = np.atleast_1d(np.asarray(mag, dtype=float))
    return {b: f0[b] * 10.0 ** (-0.4 * mag) for b in bands}


def mag0_band_fluxes(bands=BANDS):
    """Return {band: flux at AB mag 0} in photons/s/cm^2.

    Loading the LSST bandpasses and integrating the SED is relatively slow and
    ``load_lsst_bandpasses`` does not cache, so compute these normalizations
    once and reuse them (e.g. when streaming many tracts).  Flux at any
    magnitude is then ``f0[band] * 10**(-0.4*mag)``.
    """
    bps = load_lsst_bandpasses()
    sed = galsim.SED("1", wave_type="nm", flux_type="fnu")
    return {b: sed.withMagnitude(0.0, bps[b]).calculateFlux(bps[b])
            for b in bands}


# --------------------------------------------------------------------------
# 3. Write parquet files + config yaml
# --------------------------------------------------------------------------
def _galaxy_arrow_table(table, idx):
    """Build the main-file pyarrow.Table for the objects selected by `idx`."""
    f64 = pa.float64()
    f32 = pa.float32()
    cols = {
        "galaxy_id": pa.array(table["galaxy_id"][idx], type=pa.int64()),
        "ra": pa.array(table["ra"][idx], type=f64),
        "dec": pa.array(table["dec"][idx], type=f64),
        "redshift": pa.array(table["redshift"][idx], type=f64),
        "redshift_hubble": pa.array(table["redshift_hubble"][idx], type=f64),
        "shear_1": pa.array(table["shear_1"][idx], type=f64),
        "shear_2": pa.array(table["shear_2"][idx], type=f64),
        "convergence": pa.array(table["convergence"][idx], type=f64),
        "size_disk_true": pa.array(table["size_disk_true"][idx], type=f32),
        "size_minor_disk_true": pa.array(table["size_minor_disk_true"][idx], type=f32),
        "sersic_disk": pa.array(table["sersic_disk"][idx], type=f32),
        "ellipticity_1_disk_true": pa.array(table["ellipticity_1_disk_true"][idx], type=f64),
        "ellipticity_2_disk_true": pa.array(table["ellipticity_2_disk_true"][idx], type=f64),
        "sed_val_disk": pa.array([table["sed_val_disk"][i] for i in idx],
                                 type=pa.list_(f64)),
        "MW_av": pa.array(table["MW_av"][idx], type=f32),
        "MW_rv": pa.array(table["MW_rv"][idx], type=f32),
    }
    return pa.table(cols)


def _flux_arrow_table(table, fluxes, idx, bands=BANDS):
    f32 = pa.float32()
    cols = {"galaxy_id": pa.array(table["galaxy_id"][idx], type=pa.int64())}
    for b in bands:
        fb = np.asarray(fluxes[b], dtype=np.float32)
        cols[f"lsst_flux_{b}"] = pa.array(fb[idx], type=f32)
    return pa.table(cols)


def write_skycatalog(table, fluxes, out_dir, catalog_name="grid_galaxies",
                     nside=32, bands=BANDS):
    """Write galaxy_<hp>.parquet, galaxy_flux_<hp>.parquet, and skyCatalog.yaml.

    Returns the path to the top-level config yaml.
    """
    cat_dir = os.path.join(out_dir, catalog_name)
    os.makedirs(cat_dir, exist_ok=True)

    # Partition objects by healpix pixel (ring ordering, matching config).
    pix = hp.ang2pix(nside, table["ra"], table["dec"], lonlat=True, nest=False)
    for hpix in np.unique(pix):
        idx = np.where(pix == hpix)[0]
        pq.write_table(_galaxy_arrow_table(table, idx),
                       os.path.join(cat_dir, f"galaxy_{hpix}.parquet"))
        pq.write_table(_flux_arrow_table(table, fluxes, idx, bands),
                       os.path.join(cat_dir, f"galaxy_flux_{hpix}.parquet"))

    return write_skycatalog_config(cat_dir, nside=nside)


def write_skycatalog_config(cat_dir, nside=32):
    """Write skyCatalog.yaml into `cat_dir` and return its path.

    Self-contained single-dir layout: the config yaml lives alongside its
    parquet files. imSim sets skycatalog_root = dirname(config yaml), so use
    catalog_dir="." to keep root/catalog_dir pointing at that same directory.
    """
    config = {
        "skycatalog_root": os.path.abspath(cat_dir),
        "catalog_dir": ".",
        "catalog_name": "skyCatalog",
        "schema_version": "1.3.0",
        "area_partition": {"type": "healpix", "ordering": "ring", "nside": nside},
        "MW_extinction_values": {
            "a_v": {"mode": "data"},
            "r_v": {"mode": "constant", "value": 3.1},
        },
        "object_types": {
            "galaxy": {
                "data_file_type": "parquet",
                "file_template": r"galaxy_(?P<healpix>\d+).parquet",
                "flux_file_template": r"galaxy_flux_(?P<healpix>\d+).parquet",
                "spatial_model": "sersic2D",
                "sed_model": "tophat",
                "MW_extinction": "F19",
                "internal_extinction": "CCM",
                "Cosmology": COSMOLOGY,
                "tophat": {
                    "bin_parameters": ["start", "width"],
                    "bins": TOPHAT_BINS,
                    "units": "angstrom",
                },
            },
        },
    }
    config_path = os.path.join(cat_dir, "skyCatalog.yaml")
    with open(config_path, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False)
    return config_path


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Distribution parameters accept a float or a scipy spec, e.g.\n"
               "  --mag norm:loc=20,scale=0.5   --hlr uniform:loc=0.2,scale=0.8")
    p.add_argument("--ra0", type=float, default=55.0, help="grid center RA [deg]")
    p.add_argument("--dec0", type=float, default=-30.0, help="grid center Dec [deg]")
    p.add_argument("--n", type=int, default=10, help="galaxies per side")
    p.add_argument("--spacing", type=float, default=20.0, help="grid spacing [arcsec]")
    p.add_argument("--mag", type=parse_param_spec, default=20.0,
                   help="AB magnitude: float or scipy dist spec")
    p.add_argument("--hlr", type=parse_param_spec, default=0.5,
                   help="half-light radius [arcsec]: float or scipy dist spec")
    p.add_argument("--ellipticity", type=parse_param_spec, default=0.0,
                   help="intrinsic |e|: float or scipy dist spec")
    p.add_argument("--pa-deg", type=parse_param_spec, default=0.0,
                   help="intrinsic PA [deg]: float or scipy dist spec")
    p.add_argument("--g1", type=parse_param_spec, default=0.0,
                   help="applied shear g1: float or scipy dist spec")
    p.add_argument("--g2", type=parse_param_spec, default=0.0,
                   help="applied shear g2: float or scipy dist spec")
    p.add_argument("--redshift", type=parse_param_spec, default=0.1,
                   help="galaxy redshift (>0): float or scipy dist spec")
    p.add_argument("--seed", type=int, default=None,
                   help="RNG seed for distribution draws")
    p.add_argument("--nside", type=int, default=32)
    p.add_argument("--out-dir", default="/sdf/data/rubin/user/brycek/grid_galaxies/catalogs")
    p.add_argument("--name", default="grid_galaxies")
    args = p.parse_args()

    table = build_grid(args.ra0, args.dec0, args.n, args.spacing,
                       mag=args.mag, hlr=args.hlr, ellipticity=args.ellipticity,
                       pa_deg=args.pa_deg, g1=args.g1, g2=args.g2,
                       redshift=args.redshift, seed=args.seed)
    fluxes = compute_lsst_fluxes(table["mag"])
    config_path = write_skycatalog(table, fluxes, args.out_dir,
                                   catalog_name=args.name, nside=args.nside)

    mag = table["mag"]
    fr = fluxes["r"]
    print(f"Wrote {len(mag)} galaxies.")
    print(f"skyCatalog config: {config_path}")
    print(f"mag        : min={mag.min():.2f} median={np.median(mag):.2f} "
          f"max={mag.max():.2f}")
    print(f"hlr [arcsec]: min={table['size_disk_true'].min():.2f} "
          f"median={np.median(table['size_disk_true']):.2f} "
          f"max={table['size_disk_true'].max():.2f}")
    print(f"r flux [ph/s/cm^2]: min={fr.min():.4g} median={np.median(fr):.4g} "
          f"max={fr.max():.4g}")


if __name__ == "__main__":
    main()

"""Build reference-catalog input CSVs from a generated grid skyCatalog.

This is adapted from jchiang87/refcat_generation's ``make_refcat_input.py``
(which processed the OpsRehearsal-3 *star* catalogs) for the **galaxy**
skyCatalogs produced by this repo's ``run_grid_skycatalog.py`` /
``dp2_footprint_skycatalog.py``.

Differences from the upstream star version:
  * object type is ``galaxy`` (not ``star``);
  * the whole catalog is processed, pixel by pixel, instead of cone-searching a
    hardcoded list of fields;
  * magnitudes are computed through the skyCatalogs SED API
    (``get_total_observer_sed().calculateMagnitude(bp)``) for every LSST band,
    which is the same photometry path imSim uses;
  * these are noiseless grid galaxies, so no error model is applied: we write a
    "truth" refcat with small fixed position/magnitude uncertainties and no
    random scatter (so there is no dependency on the ``refcat_generation``
    package).

Environment (needs skyCatalogs + the LSST bandpasses on the path):

    source ~/setup_aos_dev.sh
    setup -k -r /sdf/home/b/brycek/u/dev-repos/skyCatalogs

Usage
-----
    python make_refcat_input.py \
        --skycatalog /path/to/<catalog_name>/skyCatalog.yaml \
        --name dp2_sv320 \
        [--bands ugrizy] [--processes 32] [--outdir .]

Produces, in ``--outdir``:
  * ``initial_refcat_mags_<name>.parquet`` -- cached per-object id/ra/dec/mags.
  * ``grid_galaxies_<name>_refmags.csv``   -- the refcat input CSV consumed by
    ``convertReferenceCatalog`` (see the README for the next steps).
"""
import os
import glob
import argparse
import warnings
from collections import defaultdict
import multiprocessing

import numpy as np
import pandas as pd
import galsim

from skycatalogs import skyCatalogs
from skycatalogs.objects.base_object import ObjectList, load_lsst_bandpasses


# Small fixed "truth" uncertainties.  The butler refcat format carries error
# columns, and downstream tasks expect them to be finite and positive, but
# these grid galaxies are noiseless so we do not add random scatter.
POS_ERR_MAS = 1.0      # per-coordinate position error, milliarcseconds
MAG_ERR = 0.001        # per-band magnitude error


lsst_bps = load_lsst_bandpasses()


def _mag0_fluxes(bands):
    """Return {band: AB-mag-0 flux} through the LSST bandpasses.

    IMPORTANT: these grid catalogs use a flat tophat SED that only sets the
    SED *shape* -- the per-object brightness lives entirely in the
    ``lsst_flux_<band>`` columns (see run_grid_skycatalog.py).  So the SED
    returned by ``get_total_observer_sed`` is *not* normalized, and
    ``sed.calculateMagnitude(bp)`` would be meaningless.  Instead we read the
    stored fluxes and invert the generator's exact relation
    ``flux = f0 * 10**(-0.4*mag)`` using the same mag-0 normalization the
    generator used (mag0_band_fluxes in run_grid_skycatalog.py), so the refcat
    magnitudes round-trip to the input magnitudes.
    """
    sed = galsim.SED("1", wave_type="nm", flux_type="fnu")
    return {b: sed.withMagnitude(0.0, lsst_bps[b]).calculateFlux(lsst_bps[b])
            for b in bands}


# --------------------------------------------------------------------------
# Per-chunk magnitude computation (runs in worker processes)
# --------------------------------------------------------------------------
def compute_mags(imin, imax, bands, outfile):
    """Compute id/ra/dec and per-band AB mags for objects[imin:imax].

    ``objects`` is a module global set by ``_init_worker`` so the (unpicklable)
    skyCatalogs ObjectList is shared via fork, matching the upstream script's
    pattern.  Magnitudes are derived from the stored ``lsst_flux_<band>`` values
    (via ``get_LSST_flux``) rather than from the unnormalized tophat SED -- see
    ``_mag0_fluxes``.
    """
    global objects
    f0 = _mag0_fluxes(bands)
    data = defaultdict(list)
    for obj in objects[imin:imax]:
        data["object_id"].append(obj.id)
        data["ra"].append(obj.ra)
        data["dec"].append(obj.dec)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            for band in bands:
                flux = obj.get_LSST_flux(band)
                # AB mag; non-positive flux -> NaN (dropped downstream).
                mag = (-2.5 * np.log10(flux / f0[band])
                       if flux is not None and flux > 0 else np.nan)
                data[band].append(mag)
    df = pd.DataFrame(data)
    df.to_parquet(outfile)
    return outfile


def _init_worker(obj_list):
    global objects
    objects = obj_list


# --------------------------------------------------------------------------
# Catalog loading -- iterate every healpix pixel of the whole catalog
# --------------------------------------------------------------------------
def read_whole_skycatalog(skycatalog_yaml):
    """Return an ObjectList with every galaxy in the catalog.

    Rather than cone-searching hardcoded field centers, walk the healpix
    pixels that actually have galaxy data files and read each in full (no
    region filter -> no double counting across overlapping regions).
    """
    skycatalog_root = os.path.dirname(os.path.abspath(skycatalog_yaml))
    sky_cat = skyCatalogs.open_catalog(skycatalog_yaml,
                                       skycatalog_root=skycatalog_root)

    obj_list = ObjectList()
    for hp in sorted(sky_cat.hps_by_type("galaxy")):
        coll = sky_cat.get_object_type_by_hp(hp, "galaxy")
        if len(coll) > 0:
            obj_list.append_object_list(coll)
    return obj_list


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skycatalog", required=True,
                    help="path to the generated <catalog>/skyCatalog.yaml")
    ap.add_argument("--name", required=True,
                    help="label for output files, e.g. dp2_sv320")
    ap.add_argument("--bands", default="ugrizy",
                    help="LSST bands to write (default: ugrizy)")
    ap.add_argument("--processes", type=int, default=32,
                    help="worker processes for SED/magnitude computation")
    ap.add_argument("--outdir", default=".",
                    help="directory for the parquet cache and output CSV")
    args = ap.parse_args()

    bands = list(args.bands)
    os.makedirs(args.outdir, exist_ok=True)

    initial_refcat_data = os.path.join(
        args.outdir, f"initial_refcat_mags_{args.name}.parquet")

    if not os.path.isfile(initial_refcat_data):
        objects = read_whole_skycatalog(args.skycatalog)
        nobj = len(objects)
        print("# skyCatalog galaxies:", nobj)
        if nobj == 0:
            raise RuntimeError(f"no galaxies found in {args.skycatalog}")

        processes = max(1, min(args.processes, nobj))
        indexes = np.linspace(0, nobj, processes + 1, dtype=int)
        chunk_files = [
            os.path.join(args.outdir, f"refcat_mags_{args.name}_"
                                      f"{imin:09d}_{imax:09d}.parquet")
            for imin, imax in zip(indexes[:-1], indexes[1:])]

        with multiprocessing.Pool(processes=processes,
                                  initializer=_init_worker,
                                  initargs=(objects,)) as pool:
            workers = [
                pool.apply_async(compute_mags, (imin, imax, bands, outfile))
                for (imin, imax), outfile in zip(
                    zip(indexes[:-1], indexes[1:]), chunk_files)]
            pool.close()
            pool.join()
            _ = [w.get() for w in workers]  # re-raise any worker exception

        dfs = [pd.read_parquet(f) for f in chunk_files if os.path.isfile(f)]
        df = pd.concat(dfs, ignore_index=True)
        df.to_parquet(initial_refcat_data)
        for f in chunk_files:                # clean up the per-chunk temporaries
            if os.path.isfile(f):
                os.remove(f)
    else:
        df = pd.read_parquet(initial_refcat_data)

    # Assemble error columns.  No random scatter -- these are truth positions
    # and magnitudes -- just small finite uncertainties the refcat format wants.
    df["ra_err"] = POS_ERR_MAS     # milliarcseconds
    df["dec_err"] = POS_ERR_MAS    # milliarcseconds
    for band in bands:
        df[f"{band}_err"] = MAG_ERR

    # Write the refcat input CSV.  Columns match what the convertReferenceCatalog
    # config (config/grid_galaxies_config.py) expects: lsst_<band> magnitudes
    # and lsst_<band>_mag_err per band.
    output_file = os.path.join(args.outdir,
                               f"grid_galaxies_{args.name}_refmags.csv")
    header = ["object_id", "ra", "dec", "ra_err", "dec_err"]
    for band in bands:
        header += [f"lsst_{band}", f"lsst_{band}_mag_err"]

    with open(output_file, "w") as fobj:
        fobj.write(",".join(header) + "\n")
        for _, row in df.iterrows():
            fields = [f"{row.object_id}", f"{row.ra}", f"{row.dec}",
                      f"{row.ra_err}", f"{row.dec_err}"]
            for band in bands:
                fields += [f"{row[band]}", f"{row[f'{band}_err']}"]
            fobj.write(",".join(fields) + "\n")

    print(f"wrote {len(df)} objects to {output_file}")


if __name__ == "__main__":
    main()

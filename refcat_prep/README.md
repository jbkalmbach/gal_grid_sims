# Building a reference catalog from a generated grid skyCatalog

Turn a galaxy skyCatalog produced by `run_grid_skycatalog.py` /
`dp2_footprint_skycatalog.py` into a reference catalog the butler can ingest
and use (e.g. for astrometric/photometric calibration of the imSim data).

Adapted from [jchiang87/refcat_generation](https://github.com/jchiang87/refcat_generation),
which built refcats from the OpsRehearsal-3 *star* catalogs.  These catalogs
are *galaxies* with a flat tophat SED whose per-object brightness lives entirely
in the `lsst_flux_<band>` columns, so magnitudes are read from those fluxes
rather than by integrating the (unnormalized) SED.  These are noiseless grid
galaxies, so we write a "truth" refcat: positions/magnitudes straight from the
catalog, with small fixed uncertainties and no random scatter.

## Files
- `make_refcat_input.py` — read a skyCatalog and write a refcat-input CSV
  (`object_id`, `ra`, `dec`, `ra_err`, `dec_err`, and `lsst_<band>` /
  `lsst_<band>_mag_err` per band).
- `config/grid_galaxies_config.py` — config for the `convertReferenceCatalog`
  CLI mapping those CSV columns to the standard refcat format.

## Environment
```bash
source ~/setup_aos_dev.sh
setup -k -r /sdf/home/b/brycek/u/dev-repos/skyCatalogs
setup -k -r /sdf/home/b/brycek/u/dev-repos/imSim
```

## 1. Write the refcat-input CSV
```bash
python make_refcat_input.py \
    --skycatalog $SCRATCH/grid_catalogs/dp2_sv320/skyCatalog.yaml \
    --name dp2_sv320 \
    --bands ugrizy --processes 32 --outdir .
```
This iterates every healpix pixel in the catalog (no hardcoded field list),
computes AB magnitudes from the stored `lsst_flux_<band>` values through the
skyCatalogs object API, and writes:
- `initial_refcat_mags_dp2_sv320.parquet` — cached id/ra/dec/mags (delete to
  force a recompute).
- `grid_galaxies_dp2_sv320_refmags.csv` — the CSV consumed in step 2.

## 2. Convert to the standard refcat format
```bash
convertReferenceCatalog \
    ./refcat_out \
    config/grid_galaxies_config.py \
    "grid_galaxies_dp2_sv320_refmags.csv"
```
(The 3rd argument is a quoted glob, so you can pass several CSVs at once.)
This writes HTM-sharded refcat FITS files under
`refcat_out/<ref_dataset_name>/`, a `master_schema.fits`, and a
`filename_to_htm.ecsv` manifest.  `<ref_dataset_name>` is set in the config
(`grid_galaxies_20251023`); bump that suffix when you regenerate so versions
don't collide.

On success `convertReferenceCatalog` prints the exact butler commands for
step 3 (with the right paths filled in).

## 3. Ingest into the butler
Using `REPO` for your butler repo and `DM-NNNNN` for your ticket/run label
(`convertReferenceCatalog` prints these with the real paths):
```bash
butler register-dataset-type REPO grid_galaxies_20251023 SimpleCatalog htm7

butler ingest-files -t direct REPO grid_galaxies_20251023 refcats/DM-NNNNN \
    ./refcat_out/filename_to_htm.ecsv

butler collection-chain REPO --mode extend refcats refcats/DM-NNNNN
```
- `register-dataset-type` registers the dataset type once per repo (the name
  must match `ref_dataset_name` in the config).
- `ingest-files -t direct` ingests the sharded FITS in place — the
  `refcat_out/` tree must stay put. Copy it somewhere permanent first if
  `/tmp` or scratch, or drop `-t direct` to copy into the repo.
- `collection-chain` adds the run to the umbrella `refcats` collection that
  calibration tasks look in.

Verify:
```bash
butler query-datasets REPO grid_galaxies_20251023 --collections refcats/DM-NNNNN
```

## Notes
- Magnitudes round-trip: a catalog generated with `--mag 20` yields refcat
  `lsst_<band>_flux` values corresponding to AB mag 20 in every band (the flat
  SED gives identical colors across bands).
- To point astrometry/photometry tasks at this refcat, set their
  `connections.*RefCat` / `ref_dataset_name` to `grid_galaxies_20251023` and
  the band mapping (`filterMap`) as needed; the flux fields are
  `lsst_<band>_flux`.
- The fixed uncertainties are set at the top of `make_refcat_input.py`
  (`POS_ERR_MAS`, `MAG_ERR`); adjust there if a task needs different values.

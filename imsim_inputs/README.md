# Running imSim at NERSC on a generated grid skyCatalog

These files run imSim on a skyCatalog produced by `run_grid_skycatalog.py` or
`dp2_footprint_skycatalog.py`, pulling visit info from a DP2 opsim database.

## Files
- `imsim-grid-run.yaml` — imSim config (galsim CLI). Based on the DP2 dev base
  config, with `input.sky_catalog.obj_types: [galaxy, star]` (the grid galaxies
  plus the sparse irregular stars our generators add for calibration) and the
  catalog path pointing at a generated `skyCatalog.yaml`.
- `opsim_visits.py` — browse the opsim `observations` table to pick a visit
  (`observationId`) by field / band / RA-Dec / night.
- `fix_seqnum.py` — post-process raw outputs to correct the `SEQNUM`/`OBSID`
  headers (see step 4).

## 1. Pick a visit from the DP2 opsim db
```bash
python opsim_visits.py --db /pscratch/sd/j/jchiang8/dp2_sims_dev/dp2_visits_fixed.db \
    --target Rubin_SV_320_-15 --band r --limit 20
# or a cone search:
python opsim_visits.py --db .../dp2_visits_fixed.db --ra 320 --dec -15 --radius 3 --band r
```
Note the `observationId` you want.

## 2. Generate the skyCatalog covering that field (example: 68-tract field)
```bash
python ../dp2_footprint_skycatalog.py \
    --tracts-file ../footprints/dp2_tracts_Rubin_SV_320_-15.txt \
    --spacing 20 --mag "norm:loc=20,scale=0.5" --hlr 0.5 --seed 42 \
    --star-density 0.15 --star-mag 22 --star-exclusion-arcsec 5 --star-seed 7 \
    --out-dir $SCRATCH/grid_catalogs --name dp2_sv320
```
This writes `$SCRATCH/grid_catalogs/dp2_sv320/skyCatalog.yaml` (+ parquet for
both the `galaxy` and `star` object types). Stars are on by default at
`--star-density 0.15` stars/arcmin² (a couple dozen per detector); pass
`--star-density 0` for a galaxy-only catalog. Star positions are uniform-random
(aperiodic, so `matchPessimisticB` finds a unique pattern) but
`--star-exclusion-arcsec` (default 5″) rejects any star landing within that
radius of a grid galaxy, keeping stars isolated. Keep that radius well below
half `--spacing` or the only room left is the hex-cell centers — a regular dual
lattice that re-breaks astrometry. The irregular stars are required for
`calibrateImage` (PSF + astrometry) to run — see `../refcat_prep/README.md`.

## 3. Run imSim (one CCD; det 94 = R22_S11)
```bash
galsim imsim-grid-run.yaml \
    input.sky_catalog.file_name=$SCRATCH/grid_catalogs/dp2_sv320/skyCatalog.yaml \
    input.opsim_data.file_name=/pscratch/sd/j/jchiang8/dp2_sims_dev/dp2_visits_fixed.db \
    input.opsim_data.visit=<observationId>
```
Output goes to `output/<visit>/` (eimage + amp files + truth catalog), per the
`output.dir` template in the YAML.

## 4. Fix the SEQNUM/OBSID headers before ingest
imSim recomputes the daily sequence number by counting exposures earlier in the
`day_obs` (`imsim/opsim_data.py`) instead of using the opsim db's own `seq_num`
column.  For the DP2 `*_fixed` dbs the stored `seq_num` already decomposes the
`observationId` exactly (`observationId = day_obs * 1e5 + seq_num`), so the
recounted value is wrong and the ingested raws get the wrong `exposure` id
(e.g. visit `2025060400350` ingests as `seq_num` 105 instead of 350).

`fix_seqnum.py` repairs already-generated raws in place, deriving the correct
value from the headers (`RUNNUM - DAYOBS * 100000`) and rewriting both `SEQNUM`
and `OBSID`:
```bash
# dry run first — report changes without writing
python fix_seqnum.py --dry-run \
    $SCRATCH/grid_galaxies/imsim_output/run_dp2_footprint_grid_mag_22/

# then fix in place
python fix_seqnum.py \
    $SCRATCH/grid_galaxies/imsim_output/run_dp2_footprint_grid_mag_22/
```
It accepts files, globs, or directories (walked recursively).  Run it **before**
ingesting the raws — the butler bakes the `exposure` id in at ingest time, so
raws already ingested with the wrong `seq_num` must be re-ingested after fixing.
Assumes single-snap visits (true for DP2).

## Notes
- The visit's band must be one for which the catalog has `lsst_flux_<band>`
  (our generators write all ugrizy for both galaxies and stars, so any band
  works).
- Stars render as point sources (`galsim.DeltaFunction`) with a flat SED; their
  total per-band flux comes from `lsst_flux_<band>` just like the galaxies, so
  the reference-catalog magnitudes round-trip exactly (see `../refcat_prep/`).
- The catalog's tracts should overlap the visit's pointing, or the CCD will be
  empty. Match the field you chose in step 1 to the tracts in step 2.
- For a full focal plane: set `input.sky_catalog.full_focal_plane: True`,
  enable `output.off_detector_photons`, and raise `output.nfiles`/`output.nproc`.
- `input.atm_psf.save_file` caches the atmospheric screens under
  `atm_psf_files/`; the first run for a visit generates them.

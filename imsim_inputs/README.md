# Running imSim at NERSC on a generated grid skyCatalog

These files run imSim on a skyCatalog produced by `run_grid_skycatalog.py` or
`dp2_footprint_skycatalog.py`, pulling visit info from a DP2 opsim database.

## Files
- `imsim-grid-run.yaml` — imSim config (galsim CLI). Based on the DP2 dev base
  config, with `input.sky_catalog.obj_types: [galaxy]` (what our generators
  write) and the catalog path pointing at a generated `skyCatalog.yaml`.
- `opsim_visits.py` — browse the opsim `observations` table to pick a visit
  (`observationId`) by field / band / RA-Dec / night.

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
    --out-dir $SCRATCH/grid_catalogs --name dp2_sv320
```
This writes `$SCRATCH/grid_catalogs/dp2_sv320/skyCatalog.yaml` (+ parquet).

## 3. Run imSim (one CCD; det 94 = R22_S11)
```bash
galsim imsim-grid-run.yaml \
    input.sky_catalog.file_name=$SCRATCH/grid_catalogs/dp2_sv320/skyCatalog.yaml \
    input.opsim_data.file_name=/pscratch/sd/j/jchiang8/dp2_sims_dev/dp2_visits_fixed.db \
    input.opsim_data.visit=<observationId>
```
Output goes to `output/<visit>/` (eimage + amp files + truth catalog), per the
`output.dir` template in the YAML.

## Notes
- The visit's band must be one for which the catalog has `lsst_flux_<band>`
  (our generators write all ugrizy, so any band works).
- The catalog's tracts should overlap the visit's pointing, or the CCD will be
  empty. Match the field you chose in step 1 to the tracts in step 2.
- For a full focal plane: set `input.sky_catalog.full_focal_plane: True`,
  enable `output.off_detector_photons`, and raise `output.nfiles`/`output.nproc`.
- `input.atm_psf.save_file` caches the atmospheric screens under
  `atm_psf_files/`; the first run for a visit generates them.

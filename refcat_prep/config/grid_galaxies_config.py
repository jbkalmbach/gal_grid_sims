# Config for ConvertReferenceCatalogTask (the `convertReferenceCatalog` CLI),
# adapted from jchiang87/refcat_generation's config/uw_stars_config.py for the
# grid-galaxy refcat CSVs written by ../make_refcat_input.py.
#
# The CSV columns are:
#   object_id, ra, dec, ra_err, dec_err,
#   lsst_<band>, lsst_<band>_mag_err   for band in ugrizy
#
# ra/dec are in degrees; ra_err/dec_err in milliarcseconds; mags are AB.
#
# Usage (see ../README.md for the full sequence):
#   convertReferenceCatalog <outdir> grid_galaxies_config.py "<csv_glob>"

import os.path

# All six LSST bands are written by make_refcat_input.py.
bands = "ugrizy"

# Name under which the converted refcat is registered in the butler.  Bump the
# suffix if you regenerate so old/new versions don't collide.
config.dataset_config.ref_dataset_name = "grid_galaxies_20251023"

config.n_processes = 32

# --- identity / astrometry columns -----------------------------------------
config.id_name = "object_id"
config.ra_name = "ra"
config.dec_name = "dec"
config.ra_err_name = "ra_err"
config.dec_err_name = "dec_err"
config.coord_err_unit = "milliarcsecond"

# --- photometry columns -----------------------------------------------------
# Magnitude column per band and its error-column mapping.
config.mag_column_list = [f"lsst_{b}" for b in bands]
config.mag_err_column_map = {f"lsst_{b}": f"lsst_{b}_mag_err" for b in bands}

# The CSV has its column-name header on the first line (line 0), which is the
# reader's default `header_lines`.  (Leave replace_missing_floats_with_nan at
# its default False: the CSV has no empty cells, and enabling it errors when the
# reader finds no masked values.)
config.file_reader.format = "csv"
config.file_reader.header_lines = 0
config.file_reader.delimiter = ","

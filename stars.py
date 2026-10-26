"""Shared helpers for a sparse, irregular population of stars (point sources).

The galaxy-grid sims in this repo lay down a *regular* hex lattice of identical
exponential galaxies.  Run through the LSST Science Pipelines, ``calibrateImage``
then fails on every detector: PSF modeling finds no point-source candidates, and
``matchPessimisticB`` locks onto a degenerate self-similar lattice solution
(~1e5 arcsec scatter).  The structural fix is to add a sparse, spatially
**irregular** population of stars spanning the field, plus a matching star
reference catalog.  The irregular pattern breaks the lattice degeneracy and gives
PSF selection real point sources to work with.

This module provides the star-specific building blocks used by both generators
(``run_grid_skycatalog.py`` and ``dp2_footprint_skycatalog.py``):

  * draw irregular (uniform-random) star positions over a sky region/box,
  * resolve a scalar/scipy/array magnitude spec (via ``run_grid_skycatalog``),
  * build the skyCatalogs ``star`` main + flux pyarrow tables, and
  * emit the ``star`` object-type config fragment for skyCatalog.yaml.

Stars are rendered by imSim as ``galsim.DeltaFunction`` point sources.  Like the
grid galaxies, their total per-band flux is carried in ``lsst_flux_<band>``
columns (``get_LSST_flux`` returns the column directly when present), so a flat
stellar SED is used only for the chromatic *shape*; the magnitudes round-trip
exactly through the reference catalog.
"""
import numpy as np
import pyarrow as pa

from run_grid_skycatalog import BANDS, resolve_param, mag0_band_fluxes

# --------------------------------------------------------------------------
# id namespace / SED
# --------------------------------------------------------------------------
# Stars share the per-tract id formula ``tract*ID_STRIDE + index`` with the
# galaxies but live in a *disjoint* numeric range so a star and a galaxy never
# emit the same ``object_id`` in an imSim truth catalog (imSim sets object_id
# from the skyCatalog object's id).  Recover a star's tract with
# ``(id - STAR_ID_NAMESPACE) // ID_STRIDE``.
ID_STRIDE = 100_000_000            # must match dp2_footprint_skycatalog.ID_STRIDE
STAR_ID_NAMESPACE = 9_000_000_000_000_000  # large int64 offset, well clear of galaxies

# Flat stellar SED (shape only; brightness comes from lsst_flux_<band>).  Path is
# relative to $SIMS_SED_LIBRARY_DIR, which StarObject._get_sed joins against.
SED_FLAT = "flatSED/sed_flat.txt.gz"

# Integer tag mixed into the star RNG seed so the star draws use a stream fully
# independent of the galaxy draws (numpy seeds must be ints, not strings).
STAR_SEED_TAG = 0xADD57A5  # "add stars"

# Default point-source surface density, stars per square arcminute.  An LSSTCam
# detector is ~0.045 deg^2 ~ 162 arcmin^2, so ~0.15/arcmin^2 gives a couple dozen
# stars per detector -- plenty for PSF/astrometry without crowding the galaxy grid.
DEFAULT_STAR_DENSITY = 0.15

# Default galaxy-exclusion radius [arcsec]: reject any candidate star within this
# distance of a galaxy node so stars stay isolated (good PSF candidates, clean
# photometry).  Keep it WELL below half the galaxy spacing (e.g. 20" grid -> use
# ~5"): if it approaches the gap radius the only room left is the hex-cell
# centers, which would rebuild a regular (dual) lattice and re-break astrometry.
DEFAULT_STAR_EXCLUSION_ARCSEC = 5.0


# --------------------------------------------------------------------------
# Irregular star positions
# --------------------------------------------------------------------------
def n_stars_for_area(area_sq_deg, density_per_sqarcmin):
    """Return the integer star count for a sky area at a given surface density.

    `area_sq_deg` is in square degrees; `density_per_sqarcmin` in stars per
    square arcminute (1 deg^2 = 3600 arcmin^2).
    """
    return int(round(area_sq_deg * 3600.0 * density_per_sqarcmin))


def _radec_to_xyz(ra, dec):
    """Unit vectors (N, 3) for ra/dec [deg]; chord distance ~= angular sep."""
    ra = np.deg2rad(np.asarray(ra, dtype=float))
    dec = np.deg2rad(np.asarray(dec, dtype=float))
    cosd = np.cos(dec)
    return np.column_stack([cosd * np.cos(ra), cosd * np.sin(ra), np.sin(dec)])


def _galaxy_tree(galaxy_ra, galaxy_dec):
    """Return a KD-tree over galaxy unit vectors, or None if no galaxies given."""
    if galaxy_ra is None or len(galaxy_ra) == 0:
        return None
    from scipy.spatial import cKDTree
    return cKDTree(_radec_to_xyz(galaxy_ra, galaxy_dec))


def _exclusion_chord(exclusion_arcsec):
    """Chord-length threshold for an angular exclusion radius [arcsec].

    Querying a unit-vector KD-tree uses Euclidean (chord) distance; the chord for
    a small angle theta is 2*sin(theta/2).  Candidates whose nearest-galaxy chord
    exceeds this are kept.
    """
    theta = np.deg2rad(exclusion_arcsec / 3600.0)
    return 2.0 * np.sin(theta / 2.0)


def _far_from_galaxies(ra, dec, tree, chord):
    """Boolean mask: True where each (ra, dec) is > `chord` from every galaxy."""
    if tree is None or chord <= 0.0:
        return np.ones(np.size(ra), dtype=bool)
    dist, _ = tree.query(_radec_to_xyz(ra, dec), k=1)
    return dist > chord


def _uniform_box(ra_c, dec_c, ra_w, dec_h, k, rng):
    """Raw uniform draw of `k` sky positions in a (ra, dec) box [deg].

    The RA extent is scaled by cos(dec) so the on-sky sampling is roughly
    uniform.  Returns (ra, dec), ra wrapped to [0, 360).
    """
    dec = dec_c + (rng.random(k) - 0.5) * dec_h
    cosd = np.cos(np.deg2rad(dec))
    ra = ra_c + (rng.random(k) - 0.5) * (ra_w / cosd)
    return ra % 360.0, dec


def _rejection_sample(n, rng, draw_batch, accept, max_oversample=100):
    """Collect `n` points, keeping only those `accept`ed.

    `draw_batch(k)` returns (ra, dec) candidate arrays; `accept(ra, dec)` returns
    a boolean mask.  Stops once `n` points are kept or ``max_oversample*n``
    candidates have been tried (so a too-aggressive filter can't loop forever).
    Returns (ra, dec) arrays of length <= n.
    """
    ras, decs = [], []
    tried = 0
    cap = max(max_oversample * n, 1)
    while len(ras) < n and tried < cap:
        k = max(n - len(ras), 1) * 2
        cand_ra, cand_dec = draw_batch(k)
        tried += cand_ra.size
        m = accept(cand_ra, cand_dec)
        for r, d in zip(cand_ra[m], cand_dec[m]):
            ras.append(r)
            decs.append(d)
            if len(ras) >= n:
                break
    return np.asarray(ras[:n]), np.asarray(decs[:n])


def draw_positions_in_box(ra_c, dec_c, ra_w, dec_h, n, rng, *,
                          galaxy_ra=None, galaxy_dec=None,
                          exclusion_arcsec=0.0):
    """Draw `n` irregular (uniform-random) sky positions within a box.

    Positions are uniform in a (ra, dec) box centered on (`ra_c`, `dec_c`) with
    full width `ra_w` and height `dec_h` [deg].  Uniform random draws are
    deliberately *not* a lattice -- the irregularity is what lets
    ``matchPessimisticB`` find a unique astrometric pattern.

    If galaxy positions and a positive `exclusion_arcsec` are given, any
    candidate within that angular radius of a galaxy is rejected (so stars stay
    isolated point-source candidates).  Keep `exclusion_arcsec` well below half
    the galaxy spacing or the only room left is the hex-cell centers -- a regular
    dual lattice that would re-break astrometry.

    Returns (ra, dec) arrays in degrees (ra wrapped to [0, 360)).
    """
    if n <= 0:
        return np.empty(0), np.empty(0)
    tree = _galaxy_tree(galaxy_ra, galaxy_dec)
    chord = _exclusion_chord(exclusion_arcsec)
    draw = lambda k: _uniform_box(ra_c, dec_c, ra_w, dec_h, k, rng)
    if tree is None:
        return draw(n)
    return _rejection_sample(
        n, rng, draw, lambda r, d: _far_from_galaxies(r, d, tree, chord))


def draw_positions_in_region(region, n, rng, *, galaxy_ra=None, galaxy_dec=None,
                             exclusion_arcsec=0.0, max_oversample=100):
    """Draw `n` irregular positions inside an sphgeom `region`.

    Rejection-samples the region's bounding box, keeping candidates that fall
    inside the region and (if galaxy positions + a positive `exclusion_arcsec`
    are given) are at least that far from every galaxy.  Used for the footprint
    generator, where each tract's non-overlapping `inner_sky_region` tiles the
    sky without double counting.
    """
    import lsst.sphgeom as sph

    if n <= 0:
        return np.empty(0), np.empty(0)
    bb = region.getBoundingBox()
    ra_c = bb.getCenter().getLon().asDegrees()
    dec_c = bb.getCenter().getLat().asDegrees()
    ra_w = bb.getWidth().asDegrees()
    dec_h = bb.getHeight().asDegrees()

    tree = _galaxy_tree(galaxy_ra, galaxy_dec)
    chord = _exclusion_chord(exclusion_arcsec)

    def accept(ra, dec):
        inside = np.fromiter(
            (region.contains(sph.UnitVector3d(sph.LonLat.fromDegrees(r, d)))
             for r, d in zip(ra, dec)), dtype=bool, count=ra.size)
        return inside & _far_from_galaxies(ra, dec, tree, chord)

    return _rejection_sample(
        n, rng, lambda k: _uniform_box(ra_c, dec_c, ra_w, dec_h, k, rng),
        accept, max_oversample=max_oversample)


# --------------------------------------------------------------------------
# skyCatalogs star tables (main + flux), built from a shared index
# --------------------------------------------------------------------------
def star_arrow_tables(ids, ra, dec, mags, f0, bands=BANDS):
    """Return (main_table, flux_table) pyarrow tables for a set of stars.

    IMPORTANT: skyCatalogs joins the main and flux files **positionally** (same
    row order, same row-group layout), not by id -- so both tables here are built
    from the same inputs in the same order.  Keep that invariant when sharding.

    Main-file columns are exactly what the skyCatalogs ``star`` reader and
    ``StarObject`` require: ``object_type`` (="star"), ``id``, ``ra``, ``dec``,
    ``magnorm``, ``sed_filepath``, and ``MW_av``/``MW_rv`` (the extinction read
    has no fallback).  ``MW_av=0`` makes extinction a no-op.  Per-band fluxes use
    the same ``f0 * 10**(-0.4*mag)`` relation as the galaxies, so the refcat
    magnitudes round-trip to the input magnitudes.
    """
    ids = np.asarray(ids, dtype=np.int64)
    ra = np.asarray(ra, dtype=np.float64)
    dec = np.asarray(dec, dtype=np.float64)
    mags = np.asarray(mags, dtype=np.float64)
    n = ids.size
    f32 = pa.float32()
    f64 = pa.float64()

    main = pa.table({
        # Homogeneous "star" column: the reader only checks the first row.
        "object_type": pa.array(["star"] * n, type=pa.string()),
        "id": pa.array(ids, type=pa.int64()),
        "ra": pa.array(ra, type=f64),
        "dec": pa.array(dec, type=f64),
        # magnorm sets the flat SED's shape normalization only; total flux comes
        # from lsst_flux_<band>.  Use the magnitude itself for a sane value.
        "magnorm": pa.array(mags, type=f64),
        "sed_filepath": pa.array([SED_FLAT] * n, type=pa.string()),
        "MW_av": pa.array(np.zeros(n, dtype=np.float32), type=f32),
        "MW_rv": pa.array(np.full(n, 3.1, dtype=np.float32), type=f32),
    })

    cols = {"id": pa.array(ids, type=pa.int64())}
    for b in bands:
        fb = (f0[b] * 10.0 ** (-0.4 * mags)).astype(np.float32)
        cols[f"lsst_flux_{b}"] = pa.array(fb, type=f32)
    flux = pa.table(cols)
    return main, flux


# --------------------------------------------------------------------------
# skyCatalog.yaml config fragment for the 'star' object type
# --------------------------------------------------------------------------
def star_config_fragment(nside=32):
    """Return the ``star`` object-type dict block for a skyCatalog.yaml.

    Mirrors the skyCatalogs star template: parquet point-source main/flux files
    named ``pointsource_<hp>.parquet`` / ``pointsource_flux_<hp>.parquet``, a
    file-based stellar SED under ``$SIMS_SED_LIBRARY_DIR``, F19 Milky-Way
    extinction, and no internal extinction (point sources).
    """
    return {
        "data_file_type": "parquet",
        "file_template": r"pointsource_(?P<healpix>\d+).parquet",
        "flux_file_template": r"pointsource_flux_(?P<healpix>\d+).parquet",
        "sed_model": "file_nm",
        "sed_file_root_env_var": "SIMS_SED_LIBRARY_DIR",
        "MW_extinction": "F19",
        "internal_extinction": "None",
        "area_partition": {"type": "healpix", "ordering": "ring", "nside": nside},
    }


__all__ = [
    "ID_STRIDE", "STAR_ID_NAMESPACE", "SED_FLAT", "STAR_SEED_TAG",
    "DEFAULT_STAR_DENSITY", "DEFAULT_STAR_EXCLUSION_ARCSEC",
    "n_stars_for_area", "draw_positions_in_box", "draw_positions_in_region",
    "star_arrow_tables", "star_config_fragment",
]

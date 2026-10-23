#!/usr/bin/env python
"""List / query visits in a DP2 opsim database, to pick a visit for imSim.

imSim's OpsimDataLoader reads visit info from the ``observations`` table of an
opsim sqlite db, keyed by ``observationId`` (see imsim/opsim_data.py).  This
helper lets you browse that table -- filter by field (target/RA-Dec), band, or
day -- and print visit ids to feed to ``input.opsim_data.visit``.

Run at NERSC where the opsim db lives, e.g.:

    python opsim_visits.py --db /pscratch/sd/j/jchiang8/dp2_sims_dev/dp2_visits_fixed.db \
        --band r --ra 320 --dec -15 --radius 3 --limit 20

Or point at a field center (RA/Dec) with a search radius, sort by seeing, etc.
"""
import argparse
import sqlite3

import numpy as np


TABLE = "observations"


def _columns(con):
    return [r[0] for r in con.execute(
        f"select name from pragma_table_info('{TABLE}')")]


def query_visits(db, band=None, ra=None, dec=None, radius=None,
                 day_obs=None, target=None, limit=50, order_by=None):
    """Return (columns, rows) from the observations table matching the filters.

    Parameters
    ----------
    db : str
        Path to the opsim sqlite database.
    band : str, optional
        Filter band (matched against the ``filter`` column).
    ra, dec, radius : float, optional
        Cone search center [deg] and radius [deg] on (fieldRA, fieldDec).
    day_obs : int, optional
        Restrict to a night (integer yyyymmdd) if the column exists.
    target : str, optional
        Substring match on a target/field-name column if one exists.
    limit : int
        Max rows to return.
    order_by : str, optional
        Column to sort ascending by (default: observationStartMJD).
    """
    con = sqlite3.connect(db)
    try:
        cols = _columns(con)
        # Choose a compact, useful projection from whatever columns exist.
        preferred = ["observationId", "fieldRA", "fieldDec", "filter",
                     "observationStartMJD", "visitExposureTime", "seeingFwhmEff",
                     "airmass", "skyBrightness", "target_name", "day_obs",
                     "science_program", "observation_reason"]
        proj = [c for c in preferred if c in cols] or cols
        where, params = [], []
        if band is not None and "filter" in cols:
            where.append("filter = ?"); params.append(band)
        if day_obs is not None and "day_obs" in cols:
            where.append("day_obs = ?"); params.append(int(day_obs))
        if target is not None:
            tcol = next((c for c in ("target_name", "science_program",
                                     "observation_reason") if c in cols), None)
            if tcol is None:
                raise ValueError("no target/field-name column in this db")
            where.append(f"{tcol} LIKE ?"); params.append(f"%{target}%")
        sql = f"select {','.join(proj)} from {TABLE}"
        if where:
            sql += " where " + " and ".join(where)
        ob = order_by or ("observationStartMJD"
                          if "observationStartMJD" in cols else proj[0])
        sql += f" order by {ob}"
        rows = [dict(zip(proj, r)) for r in con.execute(sql, params)]
    finally:
        con.close()

    # Cone search (done in python so we can use great-circle distance).
    if radius is not None and ra is not None and dec is not None:
        ra0, dec0 = np.deg2rad(ra), np.deg2rad(dec)
        kept = []
        for r in rows:
            rr, dd = np.deg2rad(r["fieldRA"]), np.deg2rad(r["fieldDec"])
            sep = np.degrees(np.arccos(np.clip(
                np.sin(dd) * np.sin(dec0)
                + np.cos(dd) * np.cos(dec0) * np.cos(rr - ra0), -1, 1)))
            if sep <= radius:
                r["sep_deg"] = round(float(sep), 3)
                kept.append(r)
        rows = kept
    return rows[:limit]


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db", required=True, help="opsim sqlite db path")
    p.add_argument("--band", help="filter band u/g/r/i/z/y")
    p.add_argument("--ra", type=float, help="cone-search center RA [deg]")
    p.add_argument("--dec", type=float, help="cone-search center Dec [deg]")
    p.add_argument("--radius", type=float, help="cone-search radius [deg]")
    p.add_argument("--day-obs", type=int, help="night as yyyymmdd")
    p.add_argument("--target", help="substring match on target/field name")
    p.add_argument("--order-by", help="column to sort by")
    p.add_argument("--limit", type=int, default=50)
    args = p.parse_args()

    try:
        rows = query_visits(args.db, band=args.band, ra=args.ra, dec=args.dec,
                            radius=args.radius, day_obs=args.day_obs,
                            target=args.target, limit=args.limit,
                            order_by=args.order_by)
    except ValueError as e:
        p.error(str(e))
    if not rows:
        print("no matching visits")
        return
    keys = list(rows[0].keys())
    print("  ".join(keys))
    for r in rows:
        print("  ".join(f"{r[k]}" for k in keys))
    print(f"\n{len(rows)} visit(s). Use an observationId as "
          f"input.opsim_data.visit.")


if __name__ == "__main__":
    main()

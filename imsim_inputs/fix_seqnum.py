#!/usr/bin/env python
"""Fix the SEQNUM / OBSID headers in imSim raw outputs.

imSim's OpsimDataLoader recomputes the daily sequence number by counting
exposures earlier in the day_obs (imsim/opsim_data.py), rather than using the
opsim db's own ``seq_num`` column.  For DP2 ``*_fixed`` dbs the stored
``seq_num`` already decomposes the ``observationId`` exactly
(``observationId = day_obs * 1e5 + seq_num``), so the recounted value is wrong
and the ingested raws get the wrong ``exposure`` id.

This script repairs already-generated raws in place by deriving the correct
sequence number from the headers themselves:

    correct_seqnum = RUNNUM - DAYOBS * 100000

where ``RUNNUM`` is the original opsim ``observationId`` (preserved verbatim by
imSim) and ``DAYOBS`` is the observing day.  It then rewrites both header cards
that encode the sequence number:

  * ``SEQNUM``  -> correct_seqnum
  * ``OBSID``   -> "<telcode>_<contrllr>_<dayobs>_<correct_seqnum:06d}"

NOTE: assumes single-snap visits (DP2).  For multi-snap visits the per-snap
sequence number is RUNNUM-derived value + snap, which is not recorded in the
primary header, so this script would need the snap index to be correct.

Usage
-----
    # dry run over a whole output tree
    python fix_seqnum.py --dry-run /path/to/imsim_output/run_.../

    # fix in place
    python fix_seqnum.py /path/to/imsim_output/run_.../

    # specific files / globs
    python fix_seqnum.py file1.fits 'dir/**/*.fits'
"""
import argparse
import glob
import os
import sys

from astropy.io import fits


def correct_seqnum_from_header(hdr):
    """Return (dayobs_int, correct_seqnum) from a header, or None if n/a."""
    if "RUNNUM" not in hdr or "DAYOBS" not in hdr:
        return None
    try:
        runnum = int(hdr["RUNNUM"])
        dayobs = int(str(hdr["DAYOBS"]))
    except (TypeError, ValueError):
        return None
    correct = runnum - dayobs * 100000
    if correct < 0:
        # observationId does not decompose as day_obs*1e5 + seq_num here.
        return None
    return dayobs, correct


def rebuild_obsid(old_obsid, dayobs, correct_seqnum):
    """Replace the trailing seqnum field of an OBSID string.

    OBSID looks like "<telcode>_<contrllr>_<dayobs>_<seqnum:06d>".
    Preserve the prefix tokens, swap in the corrected 6-digit seqnum.
    """
    if not old_obsid or "_" not in str(old_obsid):
        return f"{dayobs}_{correct_seqnum:06d}"
    parts = str(old_obsid).split("_")
    parts[-1] = f"{correct_seqnum:06d}"
    return "_".join(parts)


def fix_file(path, dry_run=True):
    """Fix SEQNUM/OBSID in every HDU of ``path`` that carries them.

    Returns a list of human-readable change descriptions (empty if no change).
    """
    changes = []
    mode = "readonly" if dry_run else "update"
    with fits.open(path, mode=mode) as hdul:
        for i, hdu in enumerate(hdul):
            hdr = hdu.header
            if "SEQNUM" not in hdr:
                continue
            res = correct_seqnum_from_header(hdr)
            if res is None:
                changes.append(f"  HDU{i}: SEQNUM present but no usable "
                               f"RUNNUM/DAYOBS -> skipped")
                continue
            dayobs, correct = res
            old_seq = int(hdr["SEQNUM"])
            if old_seq != correct:
                changes.append(f"  HDU{i}: SEQNUM {old_seq} -> {correct}")
                if not dry_run:
                    hdr["SEQNUM"] = correct
            if "OBSID" in hdr:
                old_obsid = hdr["OBSID"]
                new_obsid = rebuild_obsid(old_obsid, dayobs, correct)
                if str(old_obsid) != new_obsid:
                    changes.append(f"  HDU{i}: OBSID {old_obsid!r} "
                                   f"-> {new_obsid!r}")
                    if not dry_run:
                        hdr["OBSID"] = new_obsid
        if not dry_run and changes:
            hdul.flush()
    return changes


def iter_paths(args_paths):
    """Expand args into FITS file paths (dirs are walked recursively)."""
    for p in args_paths:
        if os.path.isdir(p):
            for root, _dirs, files in os.walk(p):
                for f in files:
                    if f.endswith((".fits", ".fits.fz", ".fz")):
                        yield os.path.join(root, f)
        elif any(ch in p for ch in "*?[") or "**" in p:
            yield from glob.glob(p, recursive=True)
        else:
            yield p


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+",
                    help="FITS files, globs, or directories to fix")
    ap.add_argument("--dry-run", action="store_true",
                    help="report changes without writing")
    args = ap.parse_args()

    n_files = n_changed = 0
    for path in iter_paths(args.paths):
        n_files += 1
        try:
            changes = fix_file(path, dry_run=args.dry_run)
        except Exception as e:  # noqa: BLE001 - report and continue
            print(f"ERROR {path}: {e}", file=sys.stderr)
            continue
        if changes:
            n_changed += 1
            tag = "WOULD FIX" if args.dry_run else "FIXED"
            print(f"{tag} {path}")
            for c in changes:
                print(c)

    verb = "would change" if args.dry_run else "changed"
    print(f"\n{n_files} file(s) scanned, {n_changed} {verb}.")


if __name__ == "__main__":
    main()

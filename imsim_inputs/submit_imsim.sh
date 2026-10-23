#!/bin/bash
# Submit one array job over all (visit, detector) pairs.
# True global throttle via %N on the single --array.
#
# Examples:
#
#   # Central-raft minus det 94 (already done), all visits in CSV
#   ./submit_imsim_global.sh \
#     --visits-csv /pscratch/sd/b/brycek/grid_galaxies/imsim_inputs/visit_list.csv \
#     --dets 90,91,92,93,95,96,97,98
#
#   # First 10 visits, specific det list, dry run
#   ./submit_imsim_global.sh --visits-csv visits.csv --first-n 10 \
#       --dets 90-98 --exclude-dets 94 --dry-run
#
#   # Explicit visits, override memory/time via extra sbatch args
#   ./submit_imsim_global.sh --visit-ids 2025060400350,2025060400351 \
#       --dets 90-98 --exclude-dets 94 -- --time=02:00:00 --mem=12G
#
#   # Lower the global concurrent cap to 40
#   ./submit_imsim_global.sh --visits-csv visits.csv --dets 90-98 --exclude-dets 94 --throttle 40

set -euo pipefail

# ---- defaults -----------------------------------------------------------
VISITS_CSV=""
VISIT_IDS=""
FIRST_N=""
DETS=""
EXCLUDE_DETS=""
THROTTLE=80
SBATCH_SCRIPT="${SBATCH_SCRIPT:-imsim_pair_array.sbatch}"
EXTRA_SBATCH_ARGS=()
DRY_RUN=0
TOTAL_DETS=189
# -------------------------------------------------------------------------

usage() {
  grep '^#' "$0" | sed 's/^# \{0,1\}//'
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --visits-csv)    VISITS_CSV="$2"; shift 2;;
    --visit-ids)     VISIT_IDS="$2"; shift 2;;
    --first-n)       FIRST_N="$2"; shift 2;;
    --dets)          DETS="$2"; shift 2;;
    --exclude-dets)  EXCLUDE_DETS="$2"; shift 2;;
    --throttle)      THROTTLE="$2"; shift 2;;
    --sbatch-script) SBATCH_SCRIPT="$2"; shift 2;;
    --dry-run)       DRY_RUN=1; shift;;
    --)              shift; EXTRA_SBATCH_ARGS+=("$@"); break;;
    -h|--help)       usage;;
    *)               echo "unknown arg: $1" >&2; usage;;
  esac
done

# ---- helper: expand "0-93,95-188,200" into sorted,unique CSV ------------
expand_spec() {
  [[ -z "${1:-}" ]] && return 0
  python3 - "$1" <<'PY'
import sys
out = set()
for p in sys.argv[1].split(','):
    p = p.strip()
    if not p:
        continue
    if '-' in p:
        a, b = p.split('-', 1)
        out.update(range(int(a), int(b) + 1))
    else:
        out.add(int(p))
print(','.join(str(x) for x in sorted(out)))
PY
}

# ---- build the detector list --------------------------------------------
[[ -z "$DETS" ]] && DETS="0-$((TOTAL_DETS - 1))"
DETS_EXPANDED=$(expand_spec "$DETS")

if [[ -n "$EXCLUDE_DETS" ]]; then
  EX=$(expand_spec "$EXCLUDE_DETS")
  DETS_EXPANDED=$(python3 - "$DETS_EXPANDED" "$EX" <<'PY'
import sys
keep = set(int(x) for x in sys.argv[1].split(',') if x)
drop = set(int(x) for x in sys.argv[2].split(',') if x)
print(','.join(str(x) for x in sorted(keep - drop)))
PY
)
fi

N_DETS_LOCAL=$(awk -F, '{print NF}' <<< "$DETS_EXPANDED")
if (( N_DETS_LOCAL == 0 )); then
  echo "No detectors to run after include/exclude." >&2
  exit 1
fi

# ---- build the visit list -----------------------------------------------
if [[ -n "$VISIT_IDS" ]]; then
  VISITS_LINES=$(echo "$VISIT_IDS" | tr ',' '\n')
elif [[ -n "$VISITS_CSV" ]]; then
  VISITS_LINES=$(awk -F, 'NR==1 && $1 !~ /^[0-9]+$/ {next} {print $1}' "$VISITS_CSV")
  if [[ -n "$FIRST_N" ]]; then
    VISITS_LINES=$(echo "$VISITS_LINES" | head -n "$FIRST_N")
  fi
else
  echo "Must supply --visits-csv or --visit-ids" >&2
  exit 1
fi

VISITS_CSV_LIST=$(echo "$VISITS_LINES" | paste -sd, -)
N_VISITS=$(echo "$VISITS_LINES" | wc -l)
TOTAL_TASKS=$(( N_VISITS * N_DETS_LOCAL ))

echo "visits=$N_VISITS dets=$N_DETS_LOCAL total_tasks=$TOTAL_TASKS throttle=$THROTTLE"
echo "dets=$DETS_EXPANDED"

# ---- export the lists so --export=ALL (sbatch default) picks them up ----
export VISITS="$VISITS_CSV_LIST"
export DETECTORS="$DETS_EXPANDED"
export N_DETS="$N_DETS_LOCAL"

# ---- submit -------------------------------------------------------------
CMD=(sbatch
     --array="0-$((TOTAL_TASKS - 1))%${THROTTLE}"
     --job-name="imsim-grid"
     "${EXTRA_SBATCH_ARGS[@]}"
     "$SBATCH_SCRIPT")

if (( DRY_RUN )); then
  printf 'DRY: '; printf '%q ' "${CMD[@]}"; printf '\n'
  echo "DRY: (would export VISITS, DETECTORS, N_DETS=$N_DETS_LOCAL)"
else
  "${CMD[@]}"
fi

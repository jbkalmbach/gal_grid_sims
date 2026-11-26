#!/bin/bash
# Submit one array job over all (visit, detector) pairs, packed onto EXCLUSIVE
# Perlmutter CPU nodes (regular QoS) instead of shared slots.
#
# Each array element is ONE whole node that runs CHUNK pairs, PACK of them
# concurrently (in successive waves) across the node's cores.
#
# Pick how the work is split one of two ways:
#   --nodes N   ask for exactly N nodes; CHUNK is derived as ceil(pairs/N).
#   --waves W   set CHUNK = PACK*W per node; node count = ceil(pairs/CHUNK).
# --nodes wins if both are given.
#
# Examples:
#
#   # Just give me 4 nodes for this sweep
#   ./submit_imsim_node.sh --visits-csv visits.csv --dets 90-98 --nodes 4
#
#
#   # Central-raft minus det 94, all visits, 60 concurrent per node, 4 waves
#   ./submit_imsim_node.sh \
#     --visits-csv /pscratch/sd/b/brycek/grid_galaxies/imsim_inputs/visit_list.csv \
#     --dets 90,91,92,93,95,96,97,98
#
#   # Dry run: see node count, chunk size, and computed --time
#   ./submit_imsim_node.sh --visits-csv visits.csv --dets 90-98 --exclude-dets 94 \
#       --pack 60 --waves 4 --dry-run
#
#   # One node over 9 dets of one visit (NNODES=1), smoke test
#   ./submit_imsim_node.sh --visits-csv visits.csv --first-n 1 --dets 90-98
#
#   # Cap concurrent NODES to 10, override time via extra sbatch args
#   #   (regular QoS max is 48h: https://docs.nersc.gov/jobs/policy/)
#   ./submit_imsim_node.sh --visits-csv visits.csv --dets 90-98 \
#       --node-throttle 10 -- --time=10:00:00

set -euo pipefail

# ---- defaults -----------------------------------------------------------
VISITS_CSV=""
VISIT_IDS=""
FIRST_N=""
DETS=""
EXCLUDE_DETS=""
PACK=60            # pairs run concurrently per node (~8G each on a 512G node)
WAVES=4            # waves per node when --nodes not given; CHUNK = PACK * WAVES
NODES=""           # ask for exactly this many nodes; chunk is derived to fit
NODE_THROTTLE=""   # max concurrent nodes (%N on the array); empty = unlimited
SBATCH_SCRIPT="${SBATCH_SCRIPT:-imsim_pair_node.sbatch}"
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
    --pack)          PACK="$2"; shift 2;;
    --waves)         WAVES="$2"; shift 2;;
    --nodes)         NODES="$2"; shift 2;;
    --node-throttle) NODE_THROTTLE="$2"; shift 2;;
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

# ---- chunking: pairs per node, node count, wall time --------------------
if (( PACK < 1 )); then echo "--pack must be >= 1" >&2; exit 1; fi

if [[ -n "$NODES" ]]; then
  # Honor a requested node count: spread the pairs evenly over that many nodes.
  # Clamp to [1, TOTAL_TASKS] (no point asking for more nodes than pairs), then
  # derive CHUNK and the implied wave count for the time estimate.
  if (( NODES < 1 )); then echo "--nodes must be >= 1" >&2; exit 1; fi
  REQ_NODES=$NODES
  (( REQ_NODES > TOTAL_TASKS )) && REQ_NODES=$TOTAL_TASKS
  CHUNK=$(( (TOTAL_TASKS + REQ_NODES - 1) / REQ_NODES ))   # ceil
  WAVES=$(( (CHUNK + PACK - 1) / PACK ))                   # ceil, for --time
else
  if (( WAVES < 1 )); then echo "--waves must be >= 1" >&2; exit 1; fi
  CHUNK=$(( PACK * WAVES ))
fi

# Actual node count always follows CHUNK, so trailing empty nodes never spawn.
NNODES=$(( (TOTAL_TASKS + CHUNK - 1) / CHUNK ))   # ceil

# Per-wave budget ~1h50m; time = WAVES*110min + 40min buffer, capped at 48h
# (regular QoS max on Perlmutter CPU; https://docs.nersc.gov/jobs/policy/).
# Format as HH:MM:00.
MAX_MIN=2880                                                # 48h
TIME_MIN=$(( WAVES * 110 + 40 ))
if (( TIME_MIN > MAX_MIN )); then
  TIME_MIN=$MAX_MIN
  MAX_WAVES=$(( (MAX_MIN - 40) / 110 ))                     # waves that fit in 48h
  MIN_NODES=$(( (TOTAL_TASKS + PACK * MAX_WAVES - 1) / (PACK * MAX_WAVES) ))
  echo "WARNING: ${WAVES} waves/node (chunk=${CHUNK}) won't finish in the 48h" >&2
  echo "         regular-QoS wall limit (~${MAX_WAVES} waves fit). Each node will" >&2
  echo "         time out mid-chunk. Use at least --nodes ${MIN_NODES} (or fewer" >&2
  echo "         --waves) so every node finishes." >&2
fi
TIME_LIMIT=$(printf '%02d:%02d:00' $((TIME_MIN / 60)) $((TIME_MIN % 60)))

ARRAY_SPEC="0-$((NNODES - 1))"
[[ -n "$NODE_THROTTLE" ]] && ARRAY_SPEC+="%${NODE_THROTTLE}"

echo "visits=$N_VISITS dets=$N_DETS_LOCAL total_pairs=$TOTAL_TASKS"
echo "pack=$PACK waves=$WAVES chunk=$CHUNK nnodes=$NNODES time=$TIME_LIMIT"
echo "array=$ARRAY_SPEC node_throttle=${NODE_THROTTLE:-none}"
echo "dets=$DETS_EXPANDED"

# ---- export the lists so --export=ALL (sbatch default) picks them up ----
export VISITS="$VISITS_CSV_LIST"
export DETECTORS="$DETS_EXPANDED"
export N_DETS="$N_DETS_LOCAL"
export CHUNK="$CHUNK"
export PACK="$PACK"

# ---- submit -------------------------------------------------------------
# --time is placed before EXTRA_SBATCH_ARGS so a trailing `-- --time=...` wins.
CMD=(sbatch
     --array="$ARRAY_SPEC"
     --time="$TIME_LIMIT"
     --job-name="imsim-node"
     "${EXTRA_SBATCH_ARGS[@]}"
     "$SBATCH_SCRIPT")

if (( DRY_RUN )); then
  printf 'DRY: '; printf '%q ' "${CMD[@]}"; printf '\n'
  echo "DRY: (would export VISITS, DETECTORS, N_DETS=$N_DETS_LOCAL, CHUNK=$CHUNK, PACK=$PACK)"
else
  "${CMD[@]}"
fi

#!/bin/bash
#SBATCH --account=m1727
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=08:00:00
#SBATCH --job-name=imsim-node
#SBATCH --output=logs/imsim-%A_%a.out
#SBATCH --error=logs/imsim-%A_%a.err
# NOTE: --array and --time are set at submission time by submit_imsim_node.sh.
# Each array element is ONE exclusive node that runs a chunk of (visit,detector)
# pairs, PACK of them at a time, in waves, across the node's cores.

set -euo pipefail
exec 2>&1

# ---- required from submit_imsim_node.sh (via exported env vars) ---------
: "${VISITS:?VISITS must be set (comma-separated list of visit ids)}"
: "${DETECTORS:?DETECTORS must be set (comma-separated list of det numbers)}"
: "${N_DETS:?N_DETS must be set (length of DETECTORS)}"
: "${CHUNK:?CHUNK must be set (number of pairs handled per node)}"
: "${PACK:?PACK must be set (number of pairs run concurrently per node)}"

# ---- user config --------------------------------------------------------
IMSIM_YAML="/pscratch/sd/b/brycek/grid_galaxies/imsim_inputs/imsim-grid-run.yaml"
WORKDIR="${WORKDIR:-$SLURM_SUBMIT_DIR}"
WRAPPER="${WRAPPER:-$WORKDIR/run_in_desc_stack.sh}"
IMAGE="${IMAGE:-lsstdesc/stack-jupyter:prod}"
# -------------------------------------------------------------------------

cd "$WORKDIR"
mkdir -p logs

# ---- thread caps (each packed task stays single-threaded) ---------------
# On an exclusive node we pin every worker to one thread so PACK concurrent
# single-core imSim runs don't oversubscribe the node.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export DESCSTACKUSERBASE=/pscratch/sd/b/brycek/add_packages

echo "[$(date -Is)] host=$(hostname) array_idx=${SLURM_ARRAY_TASK_ID} chunk=${CHUNK} pack=${PACK}"

# ---- decode which flat indices this node owns ---------------------------
IFS=',' read -r -a VIS_ARR <<< "$VISITS"
IFS=',' read -r -a DET_ARR <<< "$DETECTORS"
TOTAL=$(( ${#VIS_ARR[@]} * N_DETS ))

A=${SLURM_ARRAY_TASK_ID}
START=$(( A * CHUNK ))
END=$(( START + CHUNK ))
(( END > TOTAL )) && END=$TOTAL

if (( START >= TOTAL )); then
  echo "START=$START >= TOTAL=$TOTAL: nothing to do for this node." >&2
  exit 0
fi

echo "[$(date -Is)] node owns flat indices [${START}, ${END}) of ${TOTAL} total"

# ---- one-off sanity check (runs once per node, not per pair) ------------
echo "---- wrapper sanity check ----"
shifter --image="${IMAGE}" "$WRAPPER" bash -c '
  command -v galsim && galsim --version
  python -c "import galsim, imsim; print(galsim.__file__, galsim.__version__); print(imsim.__file__, imsim.__version__)"
'

# ---- per-pair worker ----------------------------------------------------
# Re-parses VISITS/DETECTORS from the exported CSV strings: bash arrays are
# NOT inherited by the xargs subshells, so we decode from the scalar env here.
run_one() {
  local IDX=$1
  local -a _vis _det
  IFS=',' read -r -a _vis <<< "$VISITS"
  IFS=',' read -r -a _det <<< "$DETECTORS"

  local V_IDX=$(( IDX / N_DETS ))
  local D_IDX=$(( IDX % N_DETS ))
  local VISIT_ID=${_vis[$V_IDX]}
  local DET_ID=${_det[$D_IDX]}

  local LOG="logs/imsim-${VISIT_ID}-${DET_ID}.log"
  {
    echo "[$(date -Is)] host=$(hostname) idx=${IDX} visit=${VISIT_ID} det=${DET_ID}"
    shifter --image="${IMAGE}" "$WRAPPER" \
      galsim "${IMSIM_YAML}" \
        input.opsim_data.visit="${VISIT_ID}" \
        output.det_num.first="${DET_ID}" \
        output.det_num.nitems=1
    local rc=$?
    echo "[$(date -Is)] visit ${VISIT_ID} det ${DET_ID} galsim exit=${rc}"
    exit $rc
  } > "$LOG" 2>&1
}
export -f run_one
export VISITS DETECTORS N_DETS IMSIM_YAML WRAPPER IMAGE

# ---- run this node's chunk, PACK at a time, in waves --------------------
echo "---- real run: $((END - START)) pairs, ${PACK} concurrent ----"
rc=0
seq "$START" "$((END - 1))" \
  | xargs -P "$PACK" -I{} bash -c 'run_one "$@"' _ {} \
  || rc=$?

echo "[$(date -Is)] node done: processed indices [${START}, ${END}); xargs rc=${rc}"
echo "(any non-zero worker exits are reflected in rc; see per-pair logs/imsim-<visit>-<det>.log)"
exit $rc

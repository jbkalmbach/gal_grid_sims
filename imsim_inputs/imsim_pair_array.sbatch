#!/bin/bash
#SBATCH --account=m1727
#SBATCH --constraint=cpu
#SBATCH --qos=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:50:00
#SBATCH --job-name=imsim-pair
#SBATCH --output=logs/imsim-%A_%a.out
#SBATCH --error=logs/imsim-%A_%a.err
# NOTE: --array is set at submission time by submit_imsim.sh

set -euo pipefail
exec 2>&1

# ---- required from submit_imsim_global.sh (via exported env vars) -------
: "${VISITS:?VISITS must be set (comma-separated list of visit ids)}"
: "${DETECTORS:?DETECTORS must be set (comma-separated list of det numbers)}"
: "${N_DETS:?N_DETS must be set (length of DETECTORS)}"

# ---- user config --------------------------------------------------------
IMSIM_YAML="/pscratch/sd/b/brycek/grid_galaxies/imsim_inputs/imsim-grid-run.yaml"
WORKDIR="${WORKDIR:-$SLURM_SUBMIT_DIR}"
WRAPPER="${WRAPPER:-$WORKDIR/run_in_desc_stack.sh}"
IMAGE="${IMAGE:-lsstdesc/stack-jupyter:prod}"
# -------------------------------------------------------------------------

cd "$WORKDIR"
mkdir -p logs

# ---- thread caps (prevent oversubscription of the 1-core slice) ---------
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export DESCSTACKUSERBASE=/pscratch/sd/b/brycek/add_packages

echo "image=${SHIFTER_IMAGEREQUEST:-unset}"

# ---- decode (visit, detector) from flat array index ---------------------
IFS=',' read -r -a VIS_ARR <<< "$VISITS"
IFS=',' read -r -a DET_ARR <<< "$DETECTORS"

IDX=${SLURM_ARRAY_TASK_ID}
V_IDX=$(( IDX / N_DETS ))
D_IDX=$(( IDX % N_DETS ))

if (( V_IDX < 0 || V_IDX >= ${#VIS_ARR[@]} )); then
  echo "V_IDX=$V_IDX out of range (nvisits=${#VIS_ARR[@]})" >&2
  exit 1
fi
if (( D_IDX < 0 || D_IDX >= ${#DET_ARR[@]} )); then
  echo "D_IDX=$D_IDX out of range (ndets=${#DET_ARR[@]})" >&2
  exit 1
fi

VISIT_ID=${VIS_ARR[$V_IDX]}
DET_ID=${DET_ARR[$D_IDX]}

echo "[$(date -Is)] host=$(hostname) array_idx=${IDX} v_idx=${V_IDX} d_idx=${D_IDX} visit=${VISIT_ID} det=${DET_ID}"
echo "WRAPPER=$WRAPPER"
ls -l "$WRAPPER" || echo "WRAPPER not found on host"

echo "---- shifter sanity check ----"
shifter --image="${IMAGE}" bash -c '
  echo "inside container, uname=$(uname -a)"
  echo "PATH=$PATH"
  command -v galsim || echo "galsim NOT on PATH yet (expected before wrapper setup)"
'

echo "---- wrapper sanity check ----"
shifter --image="${IMAGE}" "$WRAPPER" bash -c '
  echo "after wrapper setup"
  echo "PATH=$PATH"
  command -v galsim && galsim --version
  python -c "import galsim, imsim; print(galsim.__file__, galsim.__version__); print(imsim.__file__, imsim.__version__)"
'

echo "---- real run ----"
set -x
shifter --image="${IMAGE}" "$WRAPPER" \
  galsim "${IMSIM_YAML}" \
    input.opsim_data.visit="${VISIT_ID}" \
    output.det_num.first="${DET_ID}" \
    output.det_num.nitems=1
rc=$?
set +x
echo "galsim exit code: $rc"
# ls -l galsim.log 2>/dev/null && echo "---- galsim.log tail ----" && tail -40 galsim.log

echo "[$(date -Is)] visit ${VISIT_ID} det ${DET_ID} done"

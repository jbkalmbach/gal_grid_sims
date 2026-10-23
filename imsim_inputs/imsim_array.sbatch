#!/bin/bash
#SBATCH --account=m1727
#SBATCH --constraint=cpu
#SBATCH --qos=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --job-name=imsim-grid
#SBATCH --output=logs/imsim-%A_%a.out
#SBATCH --error=logs/imsim-%A_%a.err
# NOTE: set --array to match the number of visits in the CSV, e.g. 0-499%40
#SBATCH --array=0-257%40

set -euo pipefail
exec 2>&1

# ---- user config ---------------------------------------------------------
VISITS_CSV="/pscratch/sd/b/brycek/grid_galaxies/imsim_inputs/visit_list.csv"      # one observationId per line (header optional)
IMSIM_YAML="/pscratch/sd/b/brycek/grid_galaxies/imsim_inputs/imsim-grid-run.yaml"
WORKDIR="${WORKDIR:-$SLURM_SUBMIT_DIR}"
WRAPPER="${WRAPPER:-$WORKDIR/run_in_desc_stack.sh}"
# -------------------------------------------------------------------------

cd "$WORKDIR"
mkdir -p logs

# ---- thread caps (prevent oversubscription of the 1-core slice) ----------
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export DESCSTACKUSERBASE=/pscratch/sd/b/brycek/add_packages

echo "image=${SHIFTER_IMAGEREQUEST:-unset}"

# ---- pick the visit ------------------------------------------------------
VISIT_ID=$(
  awk -F, 'NR==1 && $1 !~ /^[0-9]+$/ {next} {print $1}' "$VISITS_CSV" \
  | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p"
)

if [[ -z "${VISIT_ID}" ]]; then
  echo "No visit id for array index ${SLURM_ARRAY_TASK_ID} in ${VISITS_CSV}" >&2
  exit 1
fi

echo "[$(date -Is)] host=$(hostname) array_idx=${SLURM_ARRAY_TASK_ID} visit=${VISIT_ID}"
echo "WRAPPER=$WRAPPER"
ls -l "$WRAPPER" || echo "WRAPPER not found on host"

echo "---- shifter sanity check ----"
shifter --image=lsstdesc/stack-jupyter:prod bash -c '
  echo "inside container, uname=$(uname -a)"
  echo "PATH=$PATH"
  command -v galsim || echo "galsim NOT on PATH yet (expected before wrapper setup)"
'

echo "---- wrapper sanity check ----"
shifter --image=lsstdesc/stack-jupyter:prod "$WRAPPER" bash -c '
  echo "after wrapper setup"
  echo "PATH=$PATH"
  command -v galsim && galsim --version
  python -c "import galsim, imsim; print(galsim.__file__, galsim.__version__); print(imsim.__file__, imsim.__version__)"
'

echo "---- real run ----"
set -x
shifter --image=lsstdesc/stack-jupyter:prod "$WRAPPER" \
  galsim "${IMSIM_YAML}" \
    input.opsim_data.visit="${VISIT_ID}" -v 2
rc=$?
set +x
echo "galsim exit code: $rc"
# ls -l galsim.log 2>/dev/null && echo "---- galsim.log tail ----" && tail -40 galsim.log

echo "[$(date -Is)] visit ${VISIT_ID} done"
#!/bin/bash -l
#SBATCH --job-name=postprocess_mimicry
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=7000
#SBATCH --time=01:00:00
#SBATCH --output=logs/postprocess_%A/slurm_%A_%a.out
#SBATCH --error=logs/postprocess_%A/slurm_%A_%a.out

set -euo pipefail

PIPELINE_DIR="${MIMICRY_PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)}}"
# shellcheck source=mimicry_pipeline.config.sh
source "${PIPELINE_DIR}/mimicry_pipeline.config.sh"

export PYTHONPATH="${MASIF_MIMICRY_SOURCE}${PYTHONPATH:+:${PYTHONPATH}}"
export PYTHONPATH="${MASIF_DB_SOURCE}:${PYTHONPATH}"

mkdir -p "${MASIF_MIMICRY_ROOT}/logs"
cd "${MASIF_MIMICRY_ROOT}"

if [[ ! -d "${TARGET_RUN_DIR}" ]]; then
  echo "Error: target run directory not found: ${TARGET_RUN_DIR}" >&2
  exit 1
fi

if [[ ! -f "${POSTPROCESS_TARGET_PDB}" ]]; then
  echo "Error: target PDB not found: ${POSTPROCESS_TARGET_PDB}" >&2
  exit 1
fi

# shellcheck source=/dev/null
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate MaSIF

# Cap BLAS/OpenMP threads to Slurm CPUs (default OMP uses all cores and can hit
# libgomp "Thread creation failed: Resource temporarily unavailable" with cpus-per-task=1).
: "${SLURM_CPUS_PER_TASK:=1}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

python -u -m masif_mimicry.postprocess \
  --target_run_dir "${TARGET_RUN_DIR}" \
  --target_pdb "${POSTPROCESS_TARGET_PDB}" \
  --database_dir "${DATABASE_DIR}" \
  --target_preprocess_dir "${TARGET_PREPROCESS_DIR_REL}" \
  --ligand "${LIGAND}" \
  --subset "${DATABASE_SUBSET_DIR}/${SLURM_ARRAY_TASK_ID}" \
  -o "${POSTPROCESS_OUT_BASENAME}_${SLURM_ARRAY_TASK_ID}.csv"

echo "Post-processing completed: ${POSTPROCESS_OUT_BASENAME}_${SLURM_ARRAY_TASK_ID}.csv"

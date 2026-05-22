#!/bin/bash -l
#SBATCH --job-name=target_sites
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/mimicry_prepare_%j.out
#SBATCH --error=logs/mimicry_prepare_%j.out

set -euo pipefail

PIPELINE_DIR="${MIMICRY_PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)}}"
# shellcheck source=mimicry_pipeline.config.sh
source "${PIPELINE_DIR}/mimicry_pipeline.config.sh"

export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}${MASIF_DB_SOURCE}:${MASIF_SEED_SOURCE}:${MASIF_MIMICRY_SOURCE}"

mkdir -p "${MASIF_MIMICRY_ROOT}/logs"
cd "${MASIF_MIMICRY_ROOT}"

apptainer exec \
  --bind "${MASIF_ROOT}:${MASIF_ROOT}" \
  --bind /scratch:/scratch \
  --bind /work:/work \
  --env PYTHONPATH="${PYTHONPATH}" \
  "${IMAGE}" \
  python -u -m masif_mimicry.define_target_sites \
    --target_preprocess_dir "${TARGET_PREPROCESS_DIR}" \
    --target_pdb "${TARGET_PDB}" \
    --target_chain "${TARGET_CHAIN}" \
    --target_residue "${TARGET_RESIDUE}" \
    --target_atom "${TARGET_ATOM}" \
    --target_ppi_id "${TARGET_PPI_ID}" \
    --num_points "${NUM_POINTS}" \
    --downsample "${TARGET_DOWNSAMPLE}" \
    --target_run_dir "${TARGET_RUN_DIR_REL}" \
    "$@"


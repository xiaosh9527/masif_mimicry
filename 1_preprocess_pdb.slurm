#!/bin/bash -l
#SBATCH --job-name=preprocess
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=7000
#SBATCH --time=00:30:00
#SBATCH --output=logs/masif_%j.out
#SBATCH --error=logs/masif_%j.out

set -euo pipefail

# SLURM runs a copy of this script from /var/spool/slurmd/...; config lives in the submit dir.
PIPELINE_DIR="${MIMICRY_PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)}}"
# shellcheck source=mimicry_pipeline.config.sh
source "${PIPELINE_DIR}/mimicry_pipeline.config.sh"

mkdir -p "${MASIF_MIMICRY_ROOT}/logs"

cd "${MASIF_MIMICRY_ROOT}"

srun apptainer exec \
  --bind "${MASIF_ROOT}:${MASIF_ROOT}" \
  "${IMAGE}" \
  bash "${MASIF_MIMICRY_ROOT}/preprocess_pdb.sh" \
    "${PDB_INPUT}" "${PDB_CHAIN}" \
    -o "${TARGET_PREPROCESS_DIR}" \
    -l "${LIGAND}" \
    -s "${LIGAND_SDF}"

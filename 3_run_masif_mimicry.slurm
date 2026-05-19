#!/bin/bash -l
#SBATCH --job-name=search_mimicry
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=02:00:00
#SBATCH --output=logs/mimicry_%A/slurm_%A_%a.out
#SBATCH --error=logs/mimicry_%A/slurm_%A_%a.out

set -euo pipefail

PIPELINE_DIR="${MIMICRY_PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)}}"
# shellcheck source=mimicry_pipeline.config.sh
source "${PIPELINE_DIR}/mimicry_pipeline.config.sh"

export PYTHONPATH="${MASIF_MIMICRY_SOURCE}${PYTHONPATH:+:${PYTHONPATH}}"
export PYTHONPATH="${MASIF_DB_SOURCE}:${MASIF_SEED_SOURCE}:${PYTHONPATH}"

mkdir -p "${MASIF_MIMICRY_ROOT}/logs/mimicry"
cd "${MASIF_MIMICRY_ROOT}"

if [[ ! -f "${TARGET_RUN_DIR}/target_sites.json" ]]; then
  echo "Error: ${TARGET_RUN_DIR}/target_sites.json not found. Run step 2 first." >&2
  exit 1
fi

apptainer exec \
  --bind /scratch:/scratch \
  --bind /work:/work \
  --env PYTHONPATH="${PYTHONPATH}" \
  "${IMAGE}" \
  python -u -m masif_mimicry.search \
    --database_dir "${DATABASE_DIR}" \
    --target_preprocess_dir "${TARGET_PREPROCESS_DIR_REL}" \
    --target_run_dir "${TARGET_RUN_DIR_REL}" \
    --top_iface_percent 0.0 \
    --count_clashes \
    --ca_clash_threshold "${CA_CLASH_THRESHOLD}" \
    --heavy_atom_clash_threshold "${HEAVY_ATOM_CLASH_THRESHOLD}" \
    --desc_dist_cutoff "${DESC_DIST_CUTOFF}" \
    --desc_dist_score_cutoff "${DESC_DIST_SCORE_CUTOFF}" \
    --split_seed_list "${DATABASE_SUBSET_DIR}/${SLURM_ARRAY_TASK_ID}"

echo "Mimicry search completed."
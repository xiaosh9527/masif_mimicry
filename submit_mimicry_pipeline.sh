#!/bin/bash
# Submit the three-stage mimicry pipeline with SLURM dependencies:
#   1_preprocess_pdb.slurm -> 2_prepare_target_sites.slurm -> 3_run_masif_mimicry.slurm
#
# Usage:
#   ./submit_mimicry_pipeline.sh
#   ./submit_mimicry_pipeline.sh --overwrite    # passed to define_target_sites.py (step 2)
#   ./submit_mimicry_pipeline.sh --skip-preprocess   # start at step 2 (prep already done)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=mimicry_pipeline.config.sh
source "${SCRIPT_DIR}/mimicry_pipeline.config.sh"

SKIP_PREPROCESS=0
PREPARE_EXTRA_ARGS=()
SEARCH_EXTRA=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip-preprocess)
      SKIP_PREPROCESS=1
      shift
      ;;
    --overwrite)
      PREPARE_EXTRA_ARGS+=(--overwrite)
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      echo "Usage: $0 [--skip-preprocess] [--overwrite]" >&2
      exit 1
      ;;
  esac
done

mkdir -p "${MASIF_MIMICRY_ROOT}/logs/mimicry"

cd "${SCRIPT_DIR}"
export MIMICRY_PIPELINE_DIR="${SCRIPT_DIR}"
SBATCH_COMMON=(--chdir="${SCRIPT_DIR}" --export=ALL)


# Step 1: preprocess PDB
if [[ "${SKIP_PREPROCESS}" -eq 0 ]]; then
  JOB1=$(sbatch --parsable "${SBATCH_COMMON[@]}" "${SCRIPT_DIR}/1_preprocess_pdb.slurm")
  echo "Step 1 (preprocess):  job ${JOB1}"
  DEP2="afterok:${JOB1}"
else
  echo "Step 1 (preprocess):  skipped"
  DEP2=""
fi

# Step 2: prepare target sites
if [[ -n "${DEP2}" ]]; then
  JOB2=$(sbatch --parsable "${SBATCH_COMMON[@]}" --dependency="${DEP2}" \
    "${SCRIPT_DIR}/2_prepare_target_sites.slurm" \
    "${PREPARE_EXTRA_ARGS[@]}")
else
  JOB2=$(sbatch --parsable "${SBATCH_COMMON[@]}" \
    "${SCRIPT_DIR}/2_prepare_target_sites.slurm" \
    "${PREPARE_EXTRA_ARGS[@]}")
fi
echo "Step 2 (target sites): job ${JOB2}"

# Step 3: run mimicry search
JOB3=$(sbatch --parsable "${SBATCH_COMMON[@]}" --dependency="afterok:${JOB2}" \
  --array="${SEARCH_ARRAY}" \
  "${SCRIPT_DIR}/3_run_masif_mimicry.slurm" \
  "${SEARCH_EXTRA[@]}")
echo "Step 3 (search):       job ${JOB3}  (array ${SEARCH_ARRAY})"


# Print summary
cat <<EOF

Pipeline submitted.
  Config: ${SCRIPT_DIR}/mimicry_pipeline.config.sh
  Target: ${TARGET_RUN_DIR}
  Seeds:  ${DATABASE_SUBSET_DIR}/<array_id>

Monitor:
  squeue -u "\$USER"
  tail -f ${MASIF_MIMICRY_ROOT}/logs/masif_*.out
  tail -f ${MASIF_MIMICRY_ROOT}/logs/mimicry_prepare_*.out
  tail -f ${MASIF_MIMICRY_ROOT}/logs/mimicry_${JOB3}_*.out
EOF

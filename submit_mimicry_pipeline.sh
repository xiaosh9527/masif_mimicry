#!/bin/bash
# Submit the MaSIF mimicry pipeline (steps controlled in mimicry_pipeline.config.sh).
#
# Usage:
#   ./submit_mimicry_pipeline.sh
#   ./submit_mimicry_pipeline.sh --overwrite   # passed to define_target_sites.py (step 2)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=mimicry_pipeline.config.sh
source "${SCRIPT_DIR}/mimicry_pipeline.config.sh"

is_true() {
  case "${1,,}" in
    true | yes | 1) return 0 ;;
    *) return 1 ;;
  esac
}

PREPARE_EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --overwrite)
      PREPARE_EXTRA_ARGS+=(--overwrite)
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      echo "Usage: $0 [--overwrite]" >&2
      echo "Enable steps via RUN_PREPROCESS, RUN_TARGET_SITES, RUN_MIMICRY_SEARCH in mimicry_pipeline.config.sh" >&2
      exit 1
      ;;
  esac
done

RUN_PREPROCESS_B=false
RUN_TARGET_SITES_B=false
RUN_MIMICRY_SEARCH_B=false
is_true "${RUN_PREPROCESS:-false}" && RUN_PREPROCESS_B=true
is_true "${RUN_TARGET_SITES:-false}" && RUN_TARGET_SITES_B=true
is_true "${RUN_MIMICRY_SEARCH:-false}" && RUN_MIMICRY_SEARCH_B=true

if ! "${RUN_PREPROCESS_B}" && ! "${RUN_TARGET_SITES_B}" && ! "${RUN_MIMICRY_SEARCH_B}"; then
  echo "Error: enable at least one of RUN_PREPROCESS, RUN_TARGET_SITES, RUN_MIMICRY_SEARCH in mimicry_pipeline.config.sh" >&2
  exit 1
fi

mkdir -p "${MASIF_MIMICRY_ROOT}/logs/mimicry"

cd "${SCRIPT_DIR}"
export MIMICRY_PIPELINE_DIR="${SCRIPT_DIR}"
SBATCH_COMMON=(--chdir="${SCRIPT_DIR}" --export=ALL)

PREV_DEP=""
JOB1=""
JOB2=""
JOB3=""

if "${RUN_PREPROCESS_B}"; then
  JOB1=$(sbatch --parsable "${SBATCH_COMMON[@]}" "${SCRIPT_DIR}/1_preprocess_pdb.slurm")
  echo "Step 1 (preprocess):   job ${JOB1}"
  PREV_DEP="afterok:${JOB1}"
else
  echo "Step 1 (preprocess):   skipped (RUN_PREPROCESS=False)"
fi

if "${RUN_TARGET_SITES_B}"; then
  if [[ -n "${PREV_DEP}" ]]; then
    JOB2=$(sbatch --parsable "${SBATCH_COMMON[@]}" --dependency="${PREV_DEP}" \
      "${SCRIPT_DIR}/2_prepare_target_sites.slurm" \
      "${PREPARE_EXTRA_ARGS[@]}")
  else
    JOB2=$(sbatch --parsable "${SBATCH_COMMON[@]}" \
      "${SCRIPT_DIR}/2_prepare_target_sites.slurm" \
      "${PREPARE_EXTRA_ARGS[@]}")
  fi
  echo "Step 2 (target sites): job ${JOB2}"
  PREV_DEP="afterok:${JOB2}"
else
  echo "Step 2 (target sites): skipped (RUN_TARGET_SITES=False)"
fi

if "${RUN_MIMICRY_SEARCH_B}"; then
  if [[ -n "${PREV_DEP}" ]]; then
    JOB3=$(sbatch --parsable "${SBATCH_COMMON[@]}" --dependency="${PREV_DEP}" \
      --array="${SEARCH_ARRAY}" \
      "${SCRIPT_DIR}/3_run_masif_mimicry.slurm")
  else
    JOB3=$(sbatch --parsable "${SBATCH_COMMON[@]}" \
      --array="${SEARCH_ARRAY}" \
      "${SCRIPT_DIR}/3_run_masif_mimicry.slurm")
  fi
  echo "Step 3 (search):       job ${JOB3}  (array ${SEARCH_ARRAY})"
else
  echo "Step 3 (search):       skipped (RUN_MIMICRY_SEARCH=False)"
fi

cat <<EOF

Pipeline submitted.
  Config: ${SCRIPT_DIR}/mimicry_pipeline.config.sh
  RUN_PREPROCESS=${RUN_PREPROCESS}  RUN_TARGET_SITES=${RUN_TARGET_SITES}  RUN_MIMICRY_SEARCH=${RUN_MIMICRY_SEARCH}
  Target: ${TARGET_RUN_DIR}
  Seeds:  ${DATABASE_SUBSET_DIR}/<array_id>

Monitor:
  squeue -u "\$USER"
  tail -f ${MASIF_MIMICRY_ROOT}/logs/masif_*.out
  tail -f ${MASIF_MIMICRY_ROOT}/logs/mimicry_prepare_*.out
EOF

if [[ -n "${JOB3}" ]]; then
  echo "  tail -f ${MASIF_MIMICRY_ROOT}/logs/mimicry_${JOB3}_*.out"
fi

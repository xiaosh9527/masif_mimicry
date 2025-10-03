#!/bin/bash
#SBATCH --job-name=proc_trunc_workflow
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/ymeng/MaSIF_mimicry_manuscript/proc_trunc_workflow_%j.out

# -------- Define paths and parameters --------

### Data paths:
WORK_DIR="/scratch/ymeng/MaSIF_mimicry_manuscript/6H0F_C_B" # Working directory to store results for filtering and truncating MaSIF results
BINDER_DIR="${WORK_DIR}/matched_proteins"                   # Directory containing the mimicry matched protein .pdb files from MaSIF mimicry 
TARGET_PDB_PATH="${WORK_DIR}/target/6H0F_B.pdb"             # Target .pdb file

### Parameters:
LIGAND="B_Y70"   # ligand chain and name
DRY_RUN=true     # generate slurm scripts only if true; submit the job array if false
TRUNC_LENGTH=86  # maximum amino acid length of the truncated structures

### Script paths:
script_dir="./bonsai_scripts"
config_file="${script_dir}/config.yaml"                                # config file storing other absolute paths and executable commands
postprocess_wrapper_script="${script_dir}/postprocess_wrapper.py"      # wrapper script to generate SLURM script for calculating postprocessing metrics on binder structures
postprocess_script="${script_dir}/postprocess_masif_mimicry.py"        # python script to calculate postprocessing metrics on binder structures
truncate_wrapper_script="${script_dir}/EvoEF2_trunc_proc_wrapper.py"   # wrapper script to generate SLURM script for truncation and rescoring
truncate_script="${script_dir}/EvoEF2_truncate.py"                     # python script to calculate optimal truncation window using EvoEF2
proc_truncate_script="${script_dir}/proc_trunc_masif_binders.py"       # python script to re-calculate metrics on truncated structures

# Helper: wait for a Slurm job (or array) to finish and validate success
wait_for_slurm_job() {
    local job_id="$1"
    local poll_seconds=10

    if ! command -v squeue >/dev/null 2>&1; then
        echo "Error: squeue not found in PATH. Cannot monitor Slurm job ${job_id}." >&2
        exit 1
    fi

    echo "Waiting for Slurm job ${job_id} to complete (polling every ${poll_seconds}s)..."
    # Wait while job (or any array task) is present in squeue
    while squeue -h -j "${job_id}" | grep -q .; do
        sleep "${poll_seconds}"
    done

    # Once it leaves the queue, query final states via sacct if available
    if command -v sacct >/dev/null 2>&1; then
        # Use pipe format to simplify parsing
        local sacct_out
        # Include children (array tasks, .batch, .extern)
        sacct_out=$(sacct -j "${job_id}" --parsable2 --noheader --format=JobID,State,ExitCode)
        echo "Final accounting for ${job_id}:"
        echo "${sacct_out}"

        # If any subjob has a failed/cancelled/timeout state, treat as failure
        if echo "${sacct_out}" | grep -E -q "\|(FAILED|CANCELLED|TIMEOUT|NODE_FAIL|PREEMPTED)\|"; then
            echo "Error: Slurm job ${job_id} did not complete successfully." >&2
            exit 1
        fi

        # Also check non-zero exit codes
        if echo "${sacct_out}" | awk -F'|' '{print $3}' | grep -E -q '^[1-9][0-9]*:|^[0-9]+:[1-9][0-9]*'; then
            echo "Error: One or more tasks in job ${job_id} reported a non-zero ExitCode." >&2
            exit 1
        fi

        echo "Job ${job_id} completed successfully."
        return 0
    fi

    # Fallback to scontrol if sacct is not available
    if command -v scontrol >/dev/null 2>&1; then
        local state
        state=$(scontrol show job "${job_id}" 2>/dev/null | awk -F= '/JobState=/{print $2}' | awk '{print $1}')
        echo "Final state for ${job_id} (scontrol): ${state}"
        case "${state}" in
            COMPLETED)
                echo "Job ${job_id} completed successfully."
                return 0
                ;;
            *)
                echo "Error: Slurm job ${job_id} ended in state ${state}." >&2
                exit 1
                ;;
        esac
    fi

    echo "Warning: Neither sacct nor scontrol available; cannot verify success for job ${job_id}. Proceeding." >&2
    return 0
}

########################################################
### Workflow begins here ###
########################################################

# -------- 1. Calculate postprocessing metrics --------
echo ""
echo "--------------------------------"
echo "Generating SLURM script for postprocessing..."
echo "--------------------------------"
echo ""

# Define parameters for postprocessing:
POSTPROCESS_DIR="${WORK_DIR}/postprocess"
postprocess_n_jobs=5

# Generate input_pdb_paths.txt and .slurm script for postprocessing:
mkdir -p ${POSTPROCESS_DIR}
postprocess_pdb_paths="${POSTPROCESS_DIR}/input_pdb_paths.txt"

# list the absolute paths of the pdb files in the binder directory
find -L "${BINDER_DIR}" -maxdepth 1 -type f -name '*.pdb' | sort > ${postprocess_pdb_paths}
echo "PDB paths for postprocessing written to: ${postprocess_pdb_paths}"

# generate the slurm script for postprocessing
python3 ${postprocess_wrapper_script} \
    --pdb_list_txt ${postprocess_pdb_paths} \
    --target_pdb_path ${TARGET_PDB_PATH} \
    --out_dir ${POSTPROCESS_DIR} \
    --ligand "${LIGAND}" \
    --process_py_path ${postprocess_script} \
    --array_size ${postprocess_n_jobs}
echo "SLURM script for postprocessing written to: ${POSTPROCESS_DIR}/postprocess_job_array.sh"

# Submit the job array if not a dry run
if [ "$DRY_RUN" = false ]; then
    submit_out=$(sbatch ${POSTPROCESS_DIR}/postprocess_job_array.sh)
    echo "${submit_out}"
    postprocess_job_id=$(echo "${submit_out}" | awk '{print $4}')
    if [ -z "${postprocess_job_id}" ]; then
        echo "Error: Failed to capture Slurm job ID for postprocessing submission." >&2
        exit 1
    fi
    wait_for_slurm_job "${postprocess_job_id}"
fi

# combine the postprocessed scores into a single csv file
postprocess_scores_csv="${POSTPROCESS_DIR}/postprocessed_scores.csv"
# Get header from first file and combine with data from all files (skipping headers)
first_file=$(ls ${POSTPROCESS_DIR}/subsets/postprocessed_scores_*.csv | head -1)
head -1 "$first_file" > ${postprocess_scores_csv}
for file in ${POSTPROCESS_DIR}/subsets/postprocessed_scores_*.csv; do
    tail -n +2 "$file" >> ${postprocess_scores_csv}
done

### Note:
### Postprocessing script identifies the binder interface residues, which are required for truncation and rescoring. 
### [Optional] Filter hits based on postprocessing metrics now, or proceed to truncation and filter at the end.


# -------- 2. Calculate optimal truncation window and re-calculate metrics on truncated structures --------
echo ""
echo "--------------------------------"
echo "Generating SLURM script for truncation and rescoring..."
echo "--------------------------------"
echo ""

# Define parameters for truncation and rescoring:
TRUNC_DIR="${WORK_DIR}/Truncate_${TRUNC_LENGTH}"
trunc_n_jobs=5

# Generate the input_csv and the .slurm script for truncation and rescoring:
mkdir -p ${TRUNC_DIR}
mkdir -p ${TRUNC_DIR}/input
trunc_input_csv="${TRUNC_DIR}/input/input_postprocess_scores.csv"
cp ${postprocess_scores_csv} ${trunc_input_csv} # use the postprocessed scores as the input for truncation and rescoring

# generate the slurm script for truncation and rescoring
python3 ${truncate_wrapper_script} \
    --root_dir ${WORK_DIR} \
    --input_csv ${trunc_input_csv} \
    --length 86 \
    --ligand ${LIGAND} \
    --n_array ${trunc_n_jobs}

# Submit the job array if not a dry run
if [ "$DRY_RUN" = false ]; then
    submit_out=$(sbatch ${TRUNC_DIR}/EvoEF2_trunc_array.sh)
    echo "${submit_out}"
    trunc_job_id=$(echo "${submit_out}" | awk '{print $4}')
    if [ -z "${trunc_job_id}" ]; then
        echo "Error: Failed to capture Slurm job ID for truncation submission." >&2
        exit 1
    fi
    wait_for_slurm_job "${trunc_job_id}"
fi

# combine the truncated score files into a single csv file
trunc_scores_csv="${TRUNC_DIR}/truncated_scores.csv"
# Get header from first file and combine with data from all files (skipping headers)
first_file=$(ls ${TRUNC_DIR}/proc_subset/proc_subset_*.csv | head -1)
head -1 "$first_file" > ${trunc_scores_csv}
for file in ${TRUNC_DIR}/proc_subset/proc_subset_*.csv; do
    tail -n +2 "$file" >> ${trunc_scores_csv}
done

#!/bin/bash

# workflow.sh: Job orchestrator for running MaSIF_seed_mimicry postprocessing 
# and truncation workflow

# 1_postprocessing.sh: launches splits input_pdb_paths.txt into N_ARRAY subsets 
# and submits a job array to calculate postprocessing metrics on matched protein
# models.

# 2_truncate.sh: launches a job array to calculate optimal truncation window using EvoEF2
# and re-calculate metrics on truncated structures.

# -------- Define paths and parameters --------

### Parameters:
LIGAND="B_Y70"   # ligand chain and name
TRUNC_LENGTH=86  # maximum amino acid length of the truncated structures

### Data paths:
WORK_DIR="/scratch/ymeng/MaSIF_mimicry_manuscript/6H0F_C_B" # Working directory to store results for filtering and truncating MaSIF results
POSTPROCESS_DIR="${WORK_DIR}/postprocess"                   # result of postprocessing metrics (default to a subdirectory of the working directory)
TRUNC_DIR="${WORK_DIR}/Truncate_${TRUNC_LENGTH}"            # result of truncation and rescoring metrics (default to a subdirectory of the working directory)
BINDER_DIR="${WORK_DIR}/matched_proteins"                   # Directory containing the mimicry matched protein .pdb files from MaSIF mimicry 
TARGET_PDB_PATH="${WORK_DIR}/target/6H0F_B.pdb"             # Target .pdb file

### Script paths: Now defined relative to repo dir
GATHER_SCRIPT="./scripts/bonsai_scripts/slurm/gather.sh"
SPLIT_SCRIPT="./scripts/bonsai_scripts/slurm/split.sh"
POSTPROCESS_SCRIPT="./scripts/bonsai_scripts/slurm/1_postprocessing.sh"
TRUNCATE_SCRIPT="./scripts/bonsai_scripts/slurm/2_truncate.sh"

### SLURM parameters:
N_ARRAY=5    # number of array jobs for postprocessing and truncation

########################################################
### Postprocessing ###
########################################################

# -------- 1. Calculate postprocessing metrics --------
echo ""
echo "--------------------------------"
echo "Generating SLURM script for postprocessing..."
echo "--------------------------------"
echo ""

# List the absolute paths of the pdb files in the binder directory as input to postprocessing script (save as input_pdb_paths.txt)
mkdir -p ${POSTPROCESS_DIR}
postprocess_pdb_paths="${POSTPROCESS_DIR}/input_pdb_paths.txt"
find -L "${BINDER_DIR}" -maxdepth 1 -type f -name '*.pdb' | sort > ${postprocess_pdb_paths}
echo "PDB paths for postprocessing written to: ${postprocess_pdb_paths}"

# submit postprocessing job array
jobid_postprocessing=$(sbatch --array 1-${N_ARRAY} ${POSTPROCESS_SCRIPT} \
                    ${postprocess_pdb_paths} ${LIGAND} ${TARGET_PDB_PATH} ${N_ARRAY} ${POSTPROCESS_DIR} | awk '{print $4}')
echo "submitted postprocessing job array with job ID: ${jobid_postprocessing}"


# submit gather job to combine the postprocessed scores into a single csv file
jobid_gather=$(sbatch --dependency=afterok:$jobid_postprocessing ${GATHER_SCRIPT} \
                    ${POSTPROCESS_DIR}/subsets "postprocessed_scores_*.csv" ${POSTPROCESS_DIR}/postprocessed_scores.csv | awk '{print $4}')
echo "submitted gather job with job ID: ${jobid_gather}"

########################################################
### Truncation ###
########################################################

### Note:
### Postprocessing script identifies the binder interface residues, which are required for truncation and rescoring. 
### [Optional] Filter hits based on postprocessing metrics now, or proceed to truncation and filter at the end.

# -------- 2. Calculate optimal truncation window and re-calculate metrics on truncated structures --------
echo ""
echo "--------------------------------"
echo "Generating SLURM script for truncation and rescoring..."
echo "--------------------------------"
echo ""

# Generate the input_csv and the .slurm script for truncation and rescoring:
mkdir -p ${TRUNC_DIR}
mkdir -p ${TRUNC_DIR}/input
# Use the gathered postprocessing scores directly as truncation input
postprocess_scores_csv="${POSTPROCESS_DIR}/postprocessed_scores.csv"
trunc_input_csv="${postprocess_scores_csv}"


# split the input_postprocess_scores.csv into N_ARRAY subsets
jobid_split=$(sbatch --dependency=afterok:$jobid_gather ${SPLIT_SCRIPT} ${trunc_input_csv} ${N_ARRAY} ${TRUNC_DIR}/input input_subset | awk '{print $4}')
echo "submitted split job with job ID: ${jobid_split}"

# Submit this job using sbatch --array=1-$N_ARRAY 2_truncate.sh <TRUNC_DIR> <TRUNC_LENGTH> <LIGAND> <N_ARRAY>
jobid_truncate=$(sbatch --dependency=afterok:$jobid_split --array 1-${N_ARRAY} ${TRUNCATE_SCRIPT} \
                    ${TRUNC_DIR} ${TRUNC_LENGTH} ${LIGAND} ${N_ARRAY} | awk '{print $4}')
echo "submitted truncation job array with job ID: ${jobid_truncate}"

# gather the truncated scores into a single csv file
jobid_gathertruncated=$(sbatch --dependency=afterok:$jobid_truncate ${GATHER_SCRIPT} \
                    ${TRUNC_DIR}/proc_subset "proc_subset_*.csv" ${TRUNC_DIR}/truncated_scores.csv | awk '{print $4}')
echo "submitted gather truncated scores job with job ID: ${jobid_gathertruncated}"



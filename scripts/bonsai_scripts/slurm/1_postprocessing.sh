#!/bin/bash -l
#SBATCH --job-name=masif_postprocessing
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=01:00:00
#SBATCH --output /scratch/%u/logs/masif_postprocessing-%A/slurm-%A_%a.out

# 1_postprocessing.sh
# SLURM script to calculate postprocessing metrics on matched protein models.
# Submit this job using sbatch --array=1-$N_ARRAY 1_postprocessing.sh <PDB_LIST> <LIGAND> <TARGET_PDB> <N_ARRAY> <OUT_DIR>
# to launch parallel array jobs.

#-------------------------------------------------------------
# Define the path to the .txt file containing pdb paths
PDB_LIST=$1
LIGAND=$2
TARGET_PDB=$3
N_ARRAY=$4 # total number of array jobs is needed to split the pdb paths into subsets
OUT_DIR=$5

root=$(git rev-parse --show-toplevel)

# Path to the postprocessing python script (relative to repo root)
POSTPROCESS_SCRIPT="$root/scripts/bonsai_scripts/python/postprocess_masif_mimicry.py"

#-------------------------------------------------------------
# Extract pdb paths corresponding to this array task
ARRAY_ID=${SLURM_ARRAY_TASK_ID}
# This selects lines where (line_number - ARRAY_ID) is divisible by $N_ARRAY
PDB_PATHS=$(awk -v id=${ARRAY_ID} ' (NR - id) % '$N_ARRAY' == 0 ' "$PDB_LIST")

# Check if there are any pdb paths for this job
if [ -z "$PDB_PATHS" ]; then
    echo "No pdb paths assigned to array task ID ${ARRAY_ID}. Exiting."
    exit 0
fi

# Echo the pdb paths that will be postprocessed
echo "Postprocessing the following .pdb paths:"
echo "$PDB_PATHS"

# Define the output CSV file, incorporating the array ID to ensure uniqueness
mkdir -p $OUT_DIR/subsets
OUTPUT_CSV=$OUT_DIR/subsets/postprocessed_scores_${ARRAY_ID}.csv

conda activate bonsai
# Run the Python script with the selected pdb paths
python3 $POSTPROCESS_SCRIPT \
    --pdb_paths $PDB_PATHS \
    -o "$OUTPUT_CSV" \
    --ligand "$LIGAND" \
    --target_pdb_path "$TARGET_PDB"

# Optional: Print a message indicating successful completion
echo "Job array task ID $ARRAY_ID completed. Output written to $OUTPUT_CSV."
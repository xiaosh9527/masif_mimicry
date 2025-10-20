#!/bin/bash -l
#SBATCH --job-name=truncate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=01:00:00
#SBATCH --output /scratch/%u/logs/truncate-%A/slurm-%A_%a.out

# 2_truncate.sh
# SLURM script to calculate optimal truncation window using EvoEF2
# and re-calculate metrics on truncated structures.
# Submit this job using sbatch --array=1-$N_ARRAY 2_truncate.sh <TRUNC_ROOT> <TRUNC_LENGTH> <LIGAND> <N_ARRAY>
# to launch parallel array jobs.

########################################
# User-defined variables and paths
########################################
TRUNC_ROOT=$1
TRUNC_LENGTH=$2
LIGAND=$3
N_ARRAY=$4

#-------------------------------------------------------------
# Path to the truncation and processing python scripts (relative to repo root)
TRUNC_SCRIPT="./scripts/bonsai_scripts/python/EvoEF2_truncate.py"
PROC_SCRIPT="./scripts/bonsai_scripts/python/proc_trunc_masif_mimicry.py"


########################################
# Generate directories and paths
########################################

TRUNC_DIR=${TRUNC_ROOT}/truncate
TEMP_DIR=${TRUNC_ROOT}/truncate/temp
TRUNC_SUBSET_DIR=${TRUNC_ROOT}/truncate/subset

#SUBSET_DIR=${TRUNC_ROOT}/subset
PROC_OUTPUT_DIR=${TRUNC_ROOT}/proc_subset

mkdir -p "$TRUNC_DIR" "$TEMP_DIR" "$TRUNC_SUBSET_DIR" "$PROC_OUTPUT_DIR"


# Get the current array task ID.
ARRAY_ID=${SLURM_ARRAY_TASK_ID}

# Input subset CSV containing only the rows for this ARRAY_ID.
INPUT_SUBSET="${TRUNC_ROOT}/input/input_subset_${ARRAY_ID}.csv"

# Truncated subset CSV with a new header column "trunc_binder_path".
TRUNC_SUBSET="${TRUNC_SUBSET_DIR}/subset_${ARRAY_ID}.csv"

# Define input and output file paths for the processing step.
PROC_TRUNC_OUTPUT="${PROC_OUTPUT_DIR}/proc_subset_${ARRAY_ID}.csv"

########################################
# Truncation Step: Build Subset CSV with truncated paths
########################################

# Read the header from the CSV and dynamically determine the column numbers.
header=$(head -n 1 "$INPUT_SUBSET")

mp_col=$(echo "$header" | awk -F, '{for(i=1;i<=NF;i++){if($i=="matched_protein_path") print i}}')
iface_col=$(echo "$header" | awk -F, '{for(i=1;i<=NF;i++){if($i=="binder_iface_residues") print i}}')
 
if [ -z "$mp_col" ] || [ -z "$iface_col" ]; then
    echo "Error: Could not find required columns 'matched_protein_path' and/or 'binder_iface_residues'."
    exit 1
fi

# Ensure output directories exist and write header with new column
mkdir -p "$(dirname "$TRUNC_SUBSET")"
echo "${header},trunc_binder_path" > "$TRUNC_SUBSET"

# Loop over each row in the input subset CSV (skipping header).
echo "Processing truncated subset: $TRUNC_SUBSET"
echo "Writing processed output to: $PROC_TRUNC_OUTPUT"

tail -n +2 "$INPUT_SUBSET" | while IFS= read -r line; do
    # Extract required columns.
    mp=$(echo "$line" | cut -d, -f"$mp_col")
    iface=$(echo "$line" | cut -d, -f"$iface_col")
    
    echo "Processing file: $mp with iface: $iface"
    
    # Run the EvoEF2_truncate.py script and capture the truncated file path.
    trunc_path=$(python3 ${TRUNC_SCRIPT} \
        --clean \
        --input "$mp" \
        --iface "$iface" \
        --length "$TRUNC_LENGTH" \
        --temp_dir "$TEMP_DIR" \
        --out_dir "$TRUNC_DIR")
    
    echo "Truncated file: $trunc_path"
    
    # Append the original row with the new truncated path to the final subset CSV.
    echo "$line,$trunc_path" >> "$TRUNC_SUBSET"
done

echo "Truncation for Array ID ${ARRAY_ID} completed. Final subset file: $TRUNC_SUBSET"


########################################
# Processing Step: Process the Truncated Subset CSV
########################################
echo "Processing truncated subset: $TRUNC_SUBSET"
echo "Writing processed output to: $PROC_TRUNC_OUTPUT"

python3 ${PROC_SCRIPT} \
    --input "${TRUNC_SUBSET}" \
    -o "${PROC_TRUNC_OUTPUT}" \
    --ligand "$LIGAND"

echo "Processing step for Array ID ${ARRAY_ID} completed. Output file: ${PROC_TRUNC_OUTPUT}"
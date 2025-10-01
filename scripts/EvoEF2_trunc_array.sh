#!/bin/bash -l
#SBATCH --job-name=EvoEF2_proc_combined
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=00:20:00
#SBATCH --array=1-500
#SBATCH --output=/scratch/ymeng/MaSIF_surface_mimicry/6H0F_C_B/Truncate_86/logs/EvoEF2_proc-%A_%a.out

########################################
# User-defined variables and paths
########################################

# Path to the CSV file with postprocessed scores.
CSV_FILE="/scratch/ymeng/MaSIF_surface_mimicry/6H0F_C_B/Truncate_86/input/filtered_postproc_scores_6H0F_C_B.csv"

# Truncation settings.
TRUNC_LENGTH=86
LIGAND="B_Y70"

# Directories for truncation.
ROOT_DIR="/scratch/ymeng/MaSIF_surface_mimicry/6H0F_C_B"

# settings from config.yaml
TRUNC_PY_PATH="/scratch/ymeng/MaSIF/Grid_search_v3.6/4_EvoEF2_truncate.py"
PROC_PY_PATH="/scratch/ymeng/MaSIF_surface_mimicry/proc_trunc_masif_mimicry.py"
MASIF_NEOSURF_ROOT_DIR="/scratch/ymeng/MaSIF/masif-neosurf"

########################################
# Generate directories and paths
########################################

TRUNC_DIR="${ROOT_DIR}/Truncate_${TRUNC_LENGTH}"
OUT_DIR="${TRUNC_DIR}/truncate"
TEMP_DIR="${OUT_DIR}/temp"
SUBSET_DIR="${OUT_DIR}/subset"
PROC_OUTPUT_DIR="${TRUNC_DIR}/proc_subset"

########################################
# Create necessary directories
########################################

mkdir -p "$OUT_DIR" "$TEMP_DIR" "$SUBSET_DIR" "$PROC_OUTPUT_DIR"

########################################
# Truncation Step: Build Subset CSV with truncated paths
########################################

# Get the current array task ID.
ARRAY_ID=${SLURM_ARRAY_TASK_ID}

# Read the header from the CSV and dynamically determine the column numbers.
header=$(head -n 1 "$CSV_FILE")

mp_col=$(echo "$header" | awk -F, '{for(i=1;i<=NF;i++){if($i=="matched_protein_path") print i}}')
iface_col=$(echo "$header" | awk -F, '{for(i=1;i<=NF;i++){if($i=="binder_iface_residues") print i}}')
 

if [ -z "$mp_col" ] || [ -z "$iface_col" ]; then
    echo "Error: Could not find required columns 'matched_protein_path' and/or 'binder_iface_residues'."
    exit 1
fi

# Create a temporary subset CSV containing only the rows for this ARRAY_ID.
TEMP_SUBSET="${SUBSET_DIR}/subset_${ARRAY_ID}_temp.csv"
{
  echo "$header"
  # The expression ((NR-2)%500 == (ARRAY_ID-1)) distributes the rows over 500 array tasks.
  awk -F, -v id="$ARRAY_ID" 'NR>1 && ((NR-2)%500 == (id-1)) {print}' "$CSV_FILE"
} > "$TEMP_SUBSET"

# Create the final subset CSV with a new header column "trunc_binder_path".
FINAL_SUBSET="${SUBSET_DIR}/subset_${ARRAY_ID}.csv"

# If the final subset already exists, skip truncation.
if [ -f "$FINAL_SUBSET" ]; then
    echo "$FINAL_SUBSET already exists. Skipping truncation step."
    rm "$TEMP_SUBSET"
else
    echo "$header,trunc_binder_path" > "$FINAL_SUBSET"
    echo "Processing subset CSV: $TEMP_SUBSET"
    echo "Writing final subset CSV with truncated paths to: $FINAL_SUBSET"

    # Loop over each row in the filtered CSV (skipping header).
    tail -n +2 "$TEMP_SUBSET" | while IFS= read -r line; do
        # Extract required columns.
        mp=$(echo "$line" | cut -d, -f"$mp_col")
        iface=$(echo "$line" | cut -d, -f"$iface_col")
        
        echo "Processing file: $mp with iface: $iface"
        
        # Run the EvoEF2_truncate.py script and capture the truncated file path.
        trunc_path=$(python3 "$TRUNC_PY_PATH"             --clean             --input "$mp"             --iface "$iface"             --length "$TRUNC_LENGTH"             --temp_dir "$TEMP_DIR"             --out_dir "$OUT_DIR")
        
        echo "Truncated file: $trunc_path"
        
        # Append the original row with the new truncated path to the final subset CSV.
        echo "$line,$trunc_path" >> "$FINAL_SUBSET"
    done

    rm "$TEMP_SUBSET"
    echo "Truncation for Array ID ${SLURM_ARRAY_TASK_ID} completed. Final subset file: $FINAL_SUBSET"
fi

########################################
# Processing Step: Process the Truncated Subset CSV
########################################

# Define input and output file paths for the processing step.
INPUT_FILE="${FINAL_SUBSET}"
OUTPUT_FILE="${PROC_OUTPUT_DIR}/proc_subset_${ARRAY_ID}.csv"


if [ -f "$OUTPUT_FILE" ]; then
    echo "$OUTPUT_FILE already exists. Done."

else
    echo "Job Array ID: ${ARRAY_ID}"
    echo "Processing input file: ${INPUT_FILE}"
    echo "Writing processed output to: ${OUTPUT_FILE}"
    
    python3 "${PROC_PY_PATH}"     --input "${INPUT_FILE}"     -o "${OUTPUT_FILE}"     --masif_neosurf_dir "${MASIF_NEOSURF_ROOT_DIR}"     --ligand "$LIGAND"

    echo "Processing step for Array ID ${ARRAY_ID} completed. Output file: ${OUTPUT_FILE}"
fi





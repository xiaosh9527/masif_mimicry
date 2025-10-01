#!/bin/bash
#SBATCH --job-name=masif_postproc
#SBATCH --output=/scratch/ymeng/MaSIF_surface_mimicry/6H0F_C_B/Postprocess/logs/slurm-%A_%a.out
#SBATCH --error=/scratch/ymeng/MaSIF_surface_mimicry/6H0F_C_B/Postprocess/logs/slurm-%A_%a.out
#SBATCH --array=0-999
#SBATCH --time=00:20:00
#SBATCH --mem=6G

# Directory containing binder .pdb files.
BINDER_DIR="/scratch/ymeng/MaSIF_surface_mimicry/6H0F_C_B/mimicry_matched_proteins"

# Get a sorted list of all binder .pdb files using find to avoid glob expansion issues.
mapfile -t ALL_PDB_FILES < <(find "${BINDER_DIR}" -maxdepth 1 -type f -name '*.pdb' | sort)

# Total number of files.
TOTAL_FILES=${#ALL_PDB_FILES[@]}

# Calculate group (chunk) size: round up total/1000.
GROUP_SIZE=$(( (TOTAL_FILES + 1000 - 1) / 1000 ))

# Calculate start index and extract the chunk for this job array task.
START_INDEX=$(( SLURM_ARRAY_TASK_ID * GROUP_SIZE ))
if [ $START_INDEX -ge $TOTAL_FILES ]; then
    echo "No files for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

# Extract the group of file paths.
GROUP_FILES=("${ALL_PDB_FILES[@]:$START_INDEX:$GROUP_SIZE}")

# Compute end index.
END_INDEX=$(( START_INDEX + ${#GROUP_FILES[@]} - 1 ))


# Build the --pdb_paths arguments.
PDB_PATHS_ARGS=""
for file in "${GROUP_FILES[@]}"; do
    PDB_PATHS_ARGS+=" $file"
done



OUT_CSV_PATH="/scratch/ymeng/MaSIF_surface_mimicry/6H0F_C_B/Postprocess/subsets/subset_${SLURM_ARRAY_TASK_ID}.csv"
if [ -f "$OUT_CSV_PATH" ]; then
    echo "$OUT_CSV_PATH already exists. Done."

else
    echo "Processing files from index ${START_INDEX} to ${END_INDEX} out of ${TOTAL_FILES}"
    # Call the python script with the group of pdb paths.
    python3 /scratch/ymeng/MaSIF_surface_mimicry/postprocess_masif_mimicry.py \
        --pdb_paths $PDB_PATHS_ARGS \
        --target_pdb_path /scratch/ymeng/MaSIF_surface_mimicry/6H0F_C_B/Preprocess/input/6H0F_B_B_Y70_502/6H0F_B.pdb \
        --query_preprocessed_root /scratch/ymeng/MaSIF_surface_mimicry/6H0F_C_B/Preprocess \
        --out_csv_file $OUT_CSV_PATH \
        --ligand B_Y70
fi




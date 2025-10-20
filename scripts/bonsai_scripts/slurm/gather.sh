#!/bin/bash -l
#SBATCH --job-name=gather
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=00:10:00
#SBATCH --output /scratch/%u/logs/gather-%A.out

# gather.sh - helper slurm script to rbind subset .csv files into one .csv file.

# Usage: sbatch gather.sh <INPUT_DIR> <PATTERN> <OUTPUT_CSV>

INPUT_DIR=$1
PATTERN=$2
OUTPUT_CSV=$3


# Ensure globs that match nothing expand to empty (not literal patterns)
shopt -s nullglob

# Expand matching files into an array (avoid using ls and handle spaces safely)
files=( "${INPUT_DIR}"/${PATTERN} )

# Abort if no files matched the pattern
if [ ${#files[@]} -eq 0 ]; then
    echo "No files found matching: ${INPUT_DIR}/${PATTERN}"
    exit 1
fi

echo "Gathering files:"
printf '%s\n' "${files[@]}"

# Get header from first file and combine with data from all files (skipping headers)
first_file="${files[0]}"
head -n 1 -- "$first_file" > "${OUTPUT_CSV}"

for file in "${files[@]}" ; do
    tail -n +2 -- "$file" >> "${OUTPUT_CSV}"
done

echo "Gathered files saved to: ${OUTPUT_CSV}"
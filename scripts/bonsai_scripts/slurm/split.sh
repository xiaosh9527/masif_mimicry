#!/bin/bash -l
#SBATCH --job-name=split
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output /scratch/%u/logs/split-%A.out

# split.sh - split a large CSV into N equal subsets, preserving the header in each.
#
# Usage:
#   sbatch split.sh <INPUT_CSV> <N_SUBSETS> <OUTPUT_DIR> <OUTPUT_BASENAME>
#
# Example:
#   sbatch split.sh /path/to/big.csv 100 /path/to/out subset
#   -> produces /path/to/out/subset_1.csv ... subset_100.csv, each with header

set -euo pipefail

INPUT_CSV=${1:-}
N_SUBSETS=${2:-}
OUTPUT_DIR=${3:-}
OUTPUT_BASENAME=${4:-subset}

if [[ -z "${INPUT_CSV}" || -z "${N_SUBSETS}" || -z "${OUTPUT_DIR}" ]]; then
    echo "Usage: sbatch split.sh <INPUT_CSV> <N_SUBSETS> <OUTPUT_DIR> <OUTPUT_BASENAME>" >&2
    exit 1
fi

if [[ ! -f "${INPUT_CSV}" ]]; then
    echo "Input CSV not found: ${INPUT_CSV}" >&2
    exit 1
fi

if ! [[ "${N_SUBSETS}" =~ ^[0-9]+$ ]] || [[ "${N_SUBSETS}" -le 0 ]]; then
    echo "N_SUBSETS must be a positive integer; got: ${N_SUBSETS}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

echo "Splitting: ${INPUT_CSV}"
echo "Number of subsets: ${N_SUBSETS}"
echo "Output dir: ${OUTPUT_DIR}"
echo "Base name: ${OUTPUT_BASENAME}"

# Read header
header=$(head -n 1 -- "${INPUT_CSV}")

# Compute total lines and body lines (exclude header)
total_lines=$(wc -l < "${INPUT_CSV}" | awk '{print $1}')
if [[ -z "${total_lines}" ]]; then total_lines=0; fi
if [[ ${total_lines} -gt 0 ]]; then
    body_lines=$(( total_lines - 1 ))
else
    body_lines=0
fi

# Pre-create output files with header
for i in $(seq 1 "${N_SUBSETS}"); do
    out_file="${OUTPUT_DIR}/${OUTPUT_BASENAME}_${i}.csv"
    printf '%s\n' "${header}" > "${out_file}"
done

# If there is no body, we are done (header-only files)
if [[ ${body_lines} -le 0 ]]; then
    echo "Input contains only header; created ${N_SUBSETS} header-only files."
    exit 0
fi

# Distribute lines as evenly as possible using a single-pass awk:
# target subset index = floor((NR-1) * N_SUBSETS / body_lines) + 1
tail -n +2 -- "${INPUT_CSV}" \
  | awk -v N="${N_SUBSETS}" -v total="${body_lines}" -v base="${OUTPUT_DIR}/${OUTPUT_BASENAME}_" '
        {
            idx = int(((NR-1) * N) / total) + 1;
            file = base idx ".csv";
            print >> file;
        }
    '

echo "Done. Created ${N_SUBSETS} files under ${OUTPUT_DIR}."



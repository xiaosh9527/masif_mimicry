#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import yaml

"""
4_EvoEF2_trunc_proc_wrapper.py
Yanxiang Meng 13.5.2025

Wrapper script to create and submit an EvoEF2 truncation/processing job.

Usage:
python3 truncate_proc_masif_binders.py \
    --root_dir /scratch/ymeng/MaSIF_search_3.6/8VLB_A_A_3JF_301 \
    --input_csv /scratch/ymeng/MaSIF_search_3.6/8VLB_A_A_3JF_301/Postprocess/postprocessed_scores_8VLB_A_A_3JF_301.csv \
    --length 86 \
    --ligand "A_3JF"\
    --n_array 500
"""

def parse_args():
    parser = argparse.ArgumentParser(
        description="Wrapper to create and submit an EvoEF2 truncation/processing job."
    )
    parser.add_argument(
        "--root_dir", required=True,
        help="Root directory (e.g., /scratch/ymeng/MaSIF_search_3.6/8VLB_A_A_3JF_301)"
    )
    parser.add_argument(
        "--input_csv", required=True,
        help="Path to the input CSV file with postprocessed scores."
    )
    parser.add_argument(
        "--length", required=True,
        help="Truncation length (e.g., 86)"
    )
    parser.add_argument(
        "--ligand", required=True,
        help="Ligand identifier (e.g., A_3JF)"
    )
    parser.add_argument(
        "--n_array", type=int, default=500,
        help="Number of array jobs for SBATCH (default: 500)"
    )
    return parser.parse_args()

def load_config():
    # Determine the directory of the current script
    script_dir = os.path.dirname(os.path.realpath(__file__))
    config_path = os.path.join(script_dir, "config.yaml")
    try:
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
    except Exception as e:
        sys.exit(f"Error reading config file {config_path}: {e}")
    return config

def create_shell_script(args, config):
    # Build paths from the config
    python_exec = config["paths"]["python_path"]  # Get python executable from config
    trunc_py_path = os.path.join(os.path.dirname(__file__), "EvoEF2_truncate.py")
    proc_py_path = os.path.join(os.path.dirname(__file__), "proc_trunc_masif_mimicry.py")

    # Create TRUNC_DIR = {root_dir}/Truncate_{length}
    trunc_dir = os.path.join(args.root_dir, f"Truncate_{args.length}")
    # Ensure the TRUNC_DIR exists
    os.makedirs(trunc_dir, exist_ok=True)

    # The shell script will be created at:
    shell_script_path = os.path.join(trunc_dir, "EvoEF2_trunc_array.sh")

    # Prepare the shell script template.
    # NOTE: All curly braces used by the shell (for variable expansion) are doubled so that
    # Python's str.format() leaves them intact.
    shell_script_template = """#!/bin/bash -l
#SBATCH --job-name=EvoEF2_proc_combined
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=00:20:00
#SBATCH --array=1-{n_array}
#SBATCH --output={root_dir}/Truncate_{trunc_length}/logs/EvoEF2_proc-%A_%a.out

########################################
# User-defined variables and paths
########################################

# Path to the CSV file with postprocessed scores.
CSV_FILE="{input_csv}"

# Truncation settings.
TRUNC_LENGTH={trunc_length}
LIGAND="{ligand}"

# Directories for truncation.
ROOT_DIR="{root_dir}"

# settings from config.yaml
TRUNC_PY_PATH="{trunc_py_path}"
PROC_PY_PATH="{proc_py_path}"

########################################
# Generate directories and paths
########################################

TRUNC_DIR="${{ROOT_DIR}}/Truncate_${{TRUNC_LENGTH}}"
OUT_DIR="${{TRUNC_DIR}}/truncate"
TEMP_DIR="${{OUT_DIR}}/temp"
SUBSET_DIR="${{OUT_DIR}}/subset"
PROC_OUTPUT_DIR="${{TRUNC_DIR}}/proc_subset"

########################################
# Create necessary directories
########################################

mkdir -p "$OUT_DIR" "$TEMP_DIR" "$SUBSET_DIR" "$PROC_OUTPUT_DIR"

########################################
# Truncation Step: Build Subset CSV with truncated paths
########################################

# Get the current array task ID.
ARRAY_ID=${{SLURM_ARRAY_TASK_ID}}

# Read the header from the CSV and dynamically determine the column numbers.
header=$(head -n 1 "$CSV_FILE")

mp_col=$(echo "$header" | awk -F, '{{for(i=1;i<=NF;i++){{if($i=="matched_protein_path") print i}}}}')
iface_col=$(echo "$header" | awk -F, '{{for(i=1;i<=NF;i++){{if($i=="binder_iface_residues") print i}}}}')
 

if [ -z "$mp_col" ] || [ -z "$iface_col" ]; then
    echo "Error: Could not find required columns 'matched_protein_path' and/or 'binder_iface_residues'."
    exit 1
fi

# Create a temporary subset CSV containing only the rows for this ARRAY_ID.
TEMP_SUBSET="${{SUBSET_DIR}}/subset_${{ARRAY_ID}}_temp.csv"
{{
  echo "$header"
  # The expression ((NR-2)%{n_array} == (ARRAY_ID-1)) distributes rows across the array tasks.
  awk -F, -v id="$ARRAY_ID" 'NR>1 && ((NR-2)%{n_array} == (id-1)) {{print}}' "$CSV_FILE"
}} > "$TEMP_SUBSET"

# Create the final subset CSV with a new header column "trunc_binder_path".
FINAL_SUBSET="${{SUBSET_DIR}}/subset_${{ARRAY_ID}}.csv"

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
        trunc_path=$({python_exec} "$TRUNC_PY_PATH" \\
            --clean \\
            --input "$mp" \\
            --iface "$iface" \\
            --length "$TRUNC_LENGTH" \\
            --temp_dir "$TEMP_DIR" \\
            --out_dir "$OUT_DIR")
        
        echo "Truncated file: $trunc_path"
        
        # Append the original row with the new truncated path to the final subset CSV.
        echo "$line,$trunc_path" >> "$FINAL_SUBSET"
    done

    rm "$TEMP_SUBSET"
    echo "Truncation for Array ID ${{SLURM_ARRAY_TASK_ID}} completed. Final subset file: $FINAL_SUBSET"
fi

########################################
# Processing Step: Process the Truncated Subset CSV
########################################

# Define input and output file paths for the processing step.
INPUT_FILE="${{FINAL_SUBSET}}"
OUTPUT_FILE="${{PROC_OUTPUT_DIR}}/proc_subset_${{ARRAY_ID}}.csv"

echo "Job Array ID: ${{ARRAY_ID}}"
echo "Processing input file: ${{INPUT_FILE}}"
echo "Writing processed output to: ${{OUTPUT_FILE}}"

{python_exec} "${{PROC_PY_PATH}}" \\
    --input "${{INPUT_FILE}}" \\
    -o "${{OUTPUT_FILE}}" \\
    --ligand "$LIGAND"

echo "Processing step for Array ID ${{ARRAY_ID}} completed. Output file: ${{OUTPUT_FILE}}"
"""

    # Fill in the template with the appropriate values.
    script_content = shell_script_template.format(
        root_dir=args.root_dir,
        input_csv=args.input_csv,
        trunc_length=args.length,
        ligand=args.ligand,
        n_array=args.n_array,
        trunc_py_path=trunc_py_path,
        proc_py_path=proc_py_path,
        python_exec=python_exec  # substitute python exec command from config
    )

    # Write the shell script to file.
    with open(shell_script_path, "w") as f:
        f.write(script_content)
    # Make the shell script executable.
    os.chmod(shell_script_path, 0o755)

    return shell_script_path


def main():
    args = parse_args()
    config = load_config()
    shell_script_path = create_shell_script(args, config)
    print(f"Shell script created at: {shell_script_path}")


if __name__ == "__main__":
    main()
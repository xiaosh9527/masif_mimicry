#!/usr/bin/env python3
"""
masif_postprocess.py

A Python wrapper script to generate a SLURM job array submission script for processing
MaSIF search outputs using multiple PDB files.

Usage:
    python3 masif_postprocess.py \
        --pdb_list_txt /path/to/postprocess_pdb_paths.txt \
        --out_dir /path/to/Postprocess/ \
        --ligand "A_3JF" \
        --target_pdb_path "/path/to/target.pdb" \
        --process_py_path "/path/to/process_search_outputs_2.2.py" \
        --job_name masif_postprocess \
        --array_size 500 \
        --mem 6G \
        --time 00:10:00 \
        --nodes 1 \
        --ntasks 1 \
        --cpus_per_task 1 \
        --script_name postprocess_job_array.sh
"""

import argparse
import os
import sys
import textwrap

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate a SLURM job array submission script for MaSIF postprocessing.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        '--pdb_list_txt',
        required=True,
        help='Path to the .txt file containing PDB paths.'
    )
    
    parser.add_argument(
        '--out_dir',
        required=True,
        help='Output directory where the submission script and outputs will be stored.'
    )
    
    parser.add_argument(
        '--ligand',
        required=True,
        help='Ligand identifier string (e.g., "A_3JF").'
    )
    
    parser.add_argument(
        '--target_pdb_path',
        required=True,
        help='Path to the target PDB file.'
    )
    
    parser.add_argument(
        '--process_py_path',
        required=True,
        help='Path to the process_search_outputs_2.2.py Python script.'
    )
    
    parser.add_argument(
        '--job_name',
        default='masif_postprocess',
        help='Name of the SLURM job array. Default: masif_postprocess'
    )
    
    parser.add_argument(
        '--array_size',
        type=int,
        default=500,
        help='Number of array tasks. Default: 500'
    )
    
    parser.add_argument(
        '--mem',
        default='6G',
        help='Memory per job. Default: 6G'
    )
    
    parser.add_argument(
        '--time',
        default='01:00:00',
        help='Time limit for each job (HH:MM:SS). Default: 01:00:00'
    )
    
    parser.add_argument(
        '--nodes',
        type=int,
        default=1,
        help='Number of nodes per job. Default: 1'
    )
    
    parser.add_argument(
        '--ntasks',
        type=int,
        default=1,
        help='Number of tasks per job. Default: 1'
    )
    
    parser.add_argument(
        '--cpus_per_task',
        type=int,
        default=1,
        help='Number of CPUs per task. Default: 1'
    )
    
    parser.add_argument(
        '--script_name',
        default='postprocess_job_array.sh',
        help='Name of the generated SLURM submission script. Default: postprocess_job_array.sh'
    )
    
    return parser.parse_args()

def validate_arguments(args):
    """Validate the provided command-line arguments."""
    # Check if pdb_list_txt exists
    if not os.path.isfile(args.pdb_list_txt):
        sys.exit(f"Error: PDB list file '{args.pdb_list_txt}' does not exist.")
    
    # Check if process_py_path exists and is a file
    if not os.path.isfile(args.process_py_path):
        sys.exit(f"Error: Python processing script '{args.process_py_path}' does not exist.")
    
    # Check if target_pdb_path exists and is a file
    if not os.path.isfile(args.target_pdb_path):
        sys.exit(f"Error: Target PDB file '{args.target_pdb_path}' does not exist.")
    
    
    
    # Create out_dir if it doesn't exist
    if not os.path.isdir(args.out_dir):
        try:
            os.makedirs(args.out_dir)
            print(f"Created output directory: {args.out_dir}")
        except Exception as e:
            sys.exit(f"Error: Could not create output directory '{args.out_dir}'. {e}")
    
    # Create logs directory inside out_dir
    logs_dir = os.path.join(args.out_dir, 'logs')
    if not os.path.isdir(logs_dir):
        try:
            os.makedirs(logs_dir)
            print(f"Created logs directory: {logs_dir}")
        except Exception as e:
            sys.exit(f"Error: Could not create logs directory '{logs_dir}'. {e}")

def generate_sh_script(args):
    """Generate the SLURM submission .sh script based on the provided arguments."""
    # Define paths
    logs_dir = os.path.join(args.out_dir, 'logs')
    script_path = os.path.join(args.out_dir, args.script_name)
    
    # Define the output CSV path pattern
    # Using ${ARRAY_ID} to substitute the array task ID
    subset_dir = os.path.join(args.out_dir, 'subsets')
    os.makedirs(subset_dir, exist_ok=True)
    output_csv_pattern = os.path.join(subset_dir, f'postprocessed_scores_${{ARRAY_ID}}.csv')
    
    # Define the SLURM output file pattern
    slurm_output_pattern = os.path.join(logs_dir, 'slurm-%A_%a.out')
    
    # Define the submission script content
    sh_script_content = textwrap.dedent(f"""\
        #!/bin/bash -l
        #SBATCH --job-name={args.job_name}
        #SBATCH --nodes={args.nodes}
        #SBATCH --ntasks={args.ntasks}
        #SBATCH --cpus-per-task={args.cpus_per_task}
        #SBATCH --mem={args.mem}
        #SBATCH --time={args.time}
        #SBATCH --array=1-{args.array_size}
        #SBATCH --output={slurm_output_pattern}
        
        # Define the path to the .txt file containing pdb paths
        PDB_LIST="{args.pdb_list_txt}"
        PY_PATH="{args.process_py_path}"
        LIGAND="{args.ligand}"
        TARGET_PDB="{args.target_pdb_path}"
        
        # Get the current array task ID
        ARRAY_ID=${{SLURM_ARRAY_TASK_ID}}
        
        # Extract pdb paths corresponding to this array task
        # This selects lines where (line_number - ARRAY_ID) is divisible by {args.array_size}
        PDB_PATHS=$(awk -v id=$ARRAY_ID ' (NR - id) % {args.array_size} == 0 ' "$PDB_LIST")
        
        # Check if there are any pdb paths for this job
        if [ -z "$PDB_PATHS" ]; then
            echo "No pdb paths assigned to array task ID $ARRAY_ID. Exiting."
            exit 0
        fi

        # Echo the pdb paths that will be postprocessed
        echo "Postprocessing the following .pdb paths:"
        echo "$PDB_PATHS"
        
        # Define the output CSV file, incorporating the array ID to ensure uniqueness
        OUTPUT_CSV="{output_csv_pattern}"
        
        # Run the Python script with the selected pdb paths
        python3 $PY_PATH \\
            --pdb_paths $PDB_PATHS \\
            -o "$OUTPUT_CSV" \\
            --ligand "$LIGAND" \\
            --target_pdb_path "$TARGET_PDB"
        
        # Optional: Print a message indicating successful completion
        echo "Job array task ID $ARRAY_ID completed. Output written to $OUTPUT_CSV."
    """)
    
    try:
        with open(script_path, 'w') as sh_file:
            sh_file.write(sh_script_content)
        print(f"SLURM submission script generated at: {script_path}")
        
        # Make the script executable
        os.chmod(script_path, 0o755)
        print("Made the script executable.")
        
    except Exception as e:
        sys.exit(f"Error: Could not write the submission script. {e}")

def main():
    """Main function to generate the SLURM submission script."""
    args = parse_arguments()
    validate_arguments(args)
    generate_sh_script(args)
    print("\nSubmission script generation completed successfully.")
    print(f"To submit the job array, run:\n\n\tsbatch {os.path.join(args.out_dir, args.script_name)}\n")

if __name__ == "__main__":
    main()
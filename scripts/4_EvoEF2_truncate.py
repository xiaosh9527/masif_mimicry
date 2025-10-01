#!/usr/bin/env python3
import os
import sys
import argparse
import shutil
import subprocess
import numpy as np
import csv
import re

"""
# 4_EvoEF2_truncate.py
# Yanxiang Meng - 25.3.2025

## Purpose
This script truncates a protein structure (PDB) file into sliding windows of a specified length and evaluates each segment using EvoEF2 to identify the most stable segment. It leverages interface residue positions and applies additional distance criteria to ensure quality truncations.

## How It Works
- **Preprocessing:** Fixes broken chain identifiers and filters out nonstandard amino acids.
- **Residue Analysis:** Determines the residue range and processes interface residues, removing outliers.
- **Window Generation:** Creates sliding windows over the interface region with a fixed length.
- **Filtering:** Discards windows where:
  - The minimum distance from any interface residue to the window’s N- or C-terminus is below a specified cutoff.
  - The distance from the mean interface residue to the nearest window terminus is below a specified cutoff.
- **Evaluation:** Runs EvoEF2 on each valid window and selects the window with the best (lowest) score.
- **Output:** Generates a CSV summary and optionally copies the best scoring PDB to an output directory.

## How to Use
Run the script using the following command:

```bash
python3 EvoEF2_truncate.py \
    --input /path/to/input.pdb \
    --iface "101_102_103_106" \
    --length 86 \
    --temp_dir /path/to/temp_dir \
    --out_dir /path/to/out_dir \
    --min_iface_dist_terminus_cutoff 2 \
    --mean_iface_dist_terminus_cutoff 5 \
    --clean
```

- **--input:** Path to the input PDB file.
- **--iface:** Interface residue numbers (underscore-separated).
- **--length:** Length of the truncation window.
- **--temp_dir:** Directory for temporary files.
- **--out_dir:** (Optional) Directory where the best PDB file will be saved.
- **--min_iface_dist_terminus_cutoff:** Minimum distance from any interface residue to the window termini.
- **--mean_iface_dist_terminus_cutoff:** Minimum distance from the mean interface residue to the nearest terminus.
- **--clean:** Remove intermediate files after processing.
"""

# define absolute paths
fix_pdb_py_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '4_fix_broken_auth_chain.py')  # in the same directory as this script
EvoEF2_dir = "~/EvoEF2"

# List of 20 standard amino acid three-letter codes
STANDARD_AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
               "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
               "THR", "TRP", "TYR", "VAL"}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Truncate a PDB based on EvoEF2 criteria.")
    parser.add_argument("--input", required=True, help="Path to input PDB file")
    parser.add_argument("--iface", required=True, help="Interface residue numbers delimited by underscore (e.g., '101_102_103_104')")
    parser.add_argument("--length", type=int, required=True, help="Window length")
    parser.add_argument("--temp_dir", required=True, help="Path to temporary directory")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("--clean", action="store_true", help="Remove intermediate files after processing")
    parser.add_argument("--out_dir", help="Destination directory to copy the best pdb file (filename remains unchanged)")
    # New criteria flags:
    parser.add_argument("--min_iface_dist_terminus_cutoff", type=int, default=2,
                        help="Minimum distance from any interface residue to either terminus of the window")
    parser.add_argument("--mean_iface_dist_terminus_cutoff", type=int, default=5,
                        help="Minimum distance from the mean interface residue to the closest terminus of the window")
    return parser.parse_args()

def safe_makedir(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)

def run_fix_broken_pdb_chain(input_pdb, output_pdb):
    cmd = ["python3", fix_pdb_py_path, "--input", input_pdb, "--output", output_pdb]
    subprocess.run(cmd, check=True)

def filter_nonstandard_amino_acids(pdb_file):
    filtered_lines = []
    with open(pdb_file, "r") as fin:
        for line in fin:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname in STANDARD_AA:
                    filtered_lines.append(line)
            else:
                filtered_lines.append(line)
    with open(pdb_file, "w") as fout:
        fout.writelines(filtered_lines)

def get_residue_range(pdb_file):
    """
    Reads the pdb_file and returns a tuple (first_resi, last_resi)
    by scanning ATOM/HETATM lines. Assumes residue number is in columns 23-26.
    """
    resi_numbers = []
    with open(pdb_file, "r") as fin:
        for line in fin:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    resi = int(line[22:26].strip())
                    resi_numbers.append(resi)
                except ValueError:
                    continue
    if not resi_numbers:
        raise ValueError("No residue numbers found in the PDB file.")
    return min(resi_numbers), max(resi_numbers)

def parse_iface(iface_str):
    return [int(x) for x in iface_str.split("_")]

def remove_outliers(iface_array, length):
    """
    Repeatedly remove the interface residue with the largest z-score until
    the range (max - min) is less than or equal to length-1.
    Returns the filtered numpy array.
    """
    iface_array = np.array(iface_array)
    while (iface_array.max() - iface_array.min()) > (length - 1):
        if iface_array.std() == 0:
            break
        zscores = np.abs((iface_array - iface_array.mean()) / iface_array.std())
        idx = np.argmax(zscores)
        iface_array = np.delete(iface_array, idx)
    return iface_array

def truncate_pdb_by_residue(input_pdb, out_pdb, start_resi, end_resi):
    with open(input_pdb, "r") as fin, open(out_pdb, "w") as fout:
        for line in fin:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    resi = int(line[22:26].strip())
                except ValueError:
                    continue
                if start_resi <= resi <= end_resi:
                    fout.write(line)
            else:
                fout.write(line)

def run_EvoEF2(pdb_path):
    cwd = os.path.expanduser(EvoEF2_dir)
    EvoEF2_exec = os.path.join(os.path.expanduser(EvoEF2_dir), "EvoEF2")
    cmd = [EvoEF2_exec, "--command=ComputeStability", f"--pdb={pdb_path}"]
    proc = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    output = proc.stdout
    total_val = None
    for line in output.splitlines():
        if "Total" in line:
            m = re.search(r"Total\s*=\s*([\-\d\.]+)", line)
            if m:
                try:
                    total_val = float(m.group(1))
                except ValueError:
                    total_val = None
                break
    if total_val is None:
        raise ValueError(f"Could not parse EvoEF2 total from output:\n{output}")
    return total_val

def main():
    args = parse_arguments()
    input_basename = os.path.splitext(os.path.basename(args.input))[0]
    base_temp_dir = os.path.join(args.temp_dir, input_basename)
    safe_makedir(base_temp_dir)

    # Run fix_broken_pdb_chain and filter nonstandard amino acids.
    fixed_pdb = os.path.join(base_temp_dir, "input.pdb")
    run_fix_broken_pdb_chain(args.input, fixed_pdb)
    filter_nonstandard_amino_acids(fixed_pdb)
    if args.verbose:
        print(f"Fixed and filtered pdb saved to {fixed_pdb}")

    # Get the residue range from the fixed pdb.
    try:
        resi_start, resi_end = get_residue_range(fixed_pdb)
    except ValueError as e:
        print(e, file=sys.stderr)
        sys.exit(1)
    pdb_length = resi_end - resi_start + 1
    if args.verbose:
        print(f"PDB residue range: {resi_start} to {resi_end} (length {pdb_length})")

    # If the pdb is shorter than --length, run EvoEF2 once on the fixed pdb.
    if pdb_length < args.length:
        if args.verbose:
            print(f"PDB length {pdb_length} is less than the specified window length {args.length}.")
            print("Running EvoEF2 on the fixed pdb without truncation.")
        evoef2_total = run_EvoEF2(fixed_pdb)
        if args.verbose:
            print(f"EvoEF2 score for fixed pdb: {evoef2_total}")
        best_file = fixed_pdb
        if args.out_dir:
            try:
                os.makedirs(args.out_dir, exist_ok=True)
                dest_filename = f"{input_basename}_{resi_start}_{resi_end}_{int(round(evoef2_total))}.pdb"
                dest_path = os.path.join(args.out_dir, dest_filename)
                shutil.copy(fixed_pdb, dest_path)
                best_file = dest_path
                if args.verbose:
                    print(f"Copied fixed pdb file to {dest_path}")
            except Exception as e:
                if args.verbose:
                    print(f"Error copying fixed pdb file: {e}", file=sys.stderr)
        if args.clean:
            try:
                shutil.rmtree(base_temp_dir)
                if args.verbose:
                    print(f"Cleaned up intermediate files in {base_temp_dir}")
            except Exception as e:
                if args.verbose:
                    print(f"Error cleaning up intermediate files: {e}", file=sys.stderr)
        print(best_file)
        return

    # Else, continue with truncation procedure.
    trunc_dir = os.path.join(base_temp_dir, "truncate")
    os.makedirs(trunc_dir, exist_ok=True)

    # Parse and filter interface residues.
    iface_array = parse_iface(args.iface)
    iface_array_filtered = remove_outliers(iface_array, args.length)
    min_iface = int(iface_array_filtered.min())
    max_iface = int(iface_array_filtered.max())
    mean_iface = np.mean(iface_array_filtered)
    if args.verbose:
        print(f"Filtered interface residues: {iface_array_filtered}")
        print(f"Interface residue range after filtering: {min_iface} to {max_iface}")
        print(f"Mean interface residue: {mean_iface}")

    # Compute candidate truncation windows
    trunc_resi_start = list(range(max_iface - args.length + 1, min_iface + 1))
    trunc_resi_end = [start + args.length - 1 for start in trunc_resi_start]

    # New: Filter candidates by checking distance criteria for each truncation window.
    valid_starts = []
    valid_ends = []
    for start, end in zip(trunc_resi_start, trunc_resi_end):
        # For each interface residue, calculate its distance to the nearest terminus.
        distances = [min(r - start, end - r) for r in iface_array_filtered]
        min_iface_dist_terminus = min(distances)
        # Calculate the distance of the mean interface residue to the nearest terminus.
        mean_iface_dist_terminus = min(mean_iface - start, end - mean_iface)
        if args.verbose:
            print(f"Window {start}-{end}: min_iface_dist_terminus={min_iface_dist_terminus}, mean_iface_dist_terminus={mean_iface_dist_terminus}")
        if (min_iface_dist_terminus >= args.min_iface_dist_terminus_cutoff) and \
           (mean_iface_dist_terminus >= args.mean_iface_dist_terminus_cutoff):
            valid_starts.append(start)
            valid_ends.append(end)
    if args.verbose:
        print(f"Valid truncation windows after applying distance criteria: {list(zip(valid_starts, valid_ends))}")

    # Process each valid truncation candidate.
    summary = []
    for start, end in zip(valid_starts, valid_ends):
        out_filename = f"{input_basename}_{start}_{end}.pdb"
        out_pdb = os.path.join(trunc_dir, out_filename)
        truncate_pdb_by_residue(fixed_pdb, out_pdb, start, end)
        # Calculate distance metrics for this window
        distances = [min(r - start, end - r) for r in iface_array_filtered]
        min_iface_dist_terminus = min(distances)
        mean_iface_dist_terminus = min(mean_iface - start, end - mean_iface)
        try:
            evoef2_total = run_EvoEF2(out_pdb)
        except Exception as e:
            if args.verbose:
                print(f"Error running EvoEF2 for {out_pdb}: {e}", file=sys.stderr)
            continue
        summary.append({
            "pdb_path": out_pdb,
            "trunc_resi_start": start,
            "trunc_resi_end": end,
            "EvoEF2": evoef2_total,
            "min_iface_dist_terminus": int(round(min_iface_dist_terminus)),
            "mean_iface_dist_terminus": int(round(mean_iface_dist_terminus))
        })
        if args.verbose:
            print(f"Processed {out_pdb} with EvoEF2 score {evoef2_total}")

    csv_path = os.path.join(base_temp_dir, f"{input_basename}_EvoEF2_trunc.csv")
    with open(csv_path, "w", newline="") as csvfile:
        fieldnames = ["pdb_path", "trunc_resi_start", "trunc_resi_end", "EvoEF2", "min_iface_dist_terminus", "mean_iface_dist_terminus"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in summary:
            row["EvoEF2"] = int(round(row["EvoEF2"]))
            writer.writerow(row)
    if args.verbose:
        print(f"Summary CSV saved to {csv_path}")

    if not summary:
        if args.verbose:
            print("No EvoEF2 results available.", file=sys.stderr)
        sys.exit(1)
    best = min(summary, key=lambda x: x["EvoEF2"])
    best_file = best["pdb_path"]
    if args.verbose:
        print(f"Best truncation: start {best['trunc_resi_start']}, end {best['trunc_resi_end']}, EvoEF2 score {int(round(best['EvoEF2']))}")

    if args.out_dir:
        try:
            os.makedirs(args.out_dir, exist_ok=True)
            dest_filename = f"{input_basename}_{best['trunc_resi_start']}_{best['trunc_resi_end']}_{int(round(best['EvoEF2']))}.pdb"
            dest_path = os.path.join(args.out_dir, dest_filename)
            shutil.copy(best["pdb_path"], dest_path)
            best_file = dest_path
            if args.verbose:
                print(f"Copied best pdb file to {dest_path}")
        except Exception as e:
            if args.verbose:
                print(f"Error copying best pdb file: {e}", file=sys.stderr)

    if args.clean:
        try:
            shutil.rmtree(base_temp_dir)
            if args.verbose:
                print(f"Cleaned up intermediate files in {base_temp_dir}")
        except Exception as e:
            if args.verbose:
                print(f"Error cleaning up intermediate files: {e}", file=sys.stderr)

    print(best_file)

if __name__ == "__main__":
    main()
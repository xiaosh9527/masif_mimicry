#!/usr/bin/env python3
import os
import sys
import argparse
import shutil
import subprocess
import numpy as np
import csv
import re
import yaml

"""
python3 EvoEF2_truncate.py \
    --clean \
    --input /path/to/input.pdb \
    --iface "101_102_103_106" \
    --length 86 \
    --temp_dir /scratch/ymeng/MaSIF_search_3.6/8VLB_A_A_3JF_301/Postprocess/temp \
    --out_dir MaSIF_search_3.6/8VLB_A_A_3JF_301/Postprocess/truncate
"""
# ------------------ Load configuration ------------------
# Adjust CONFIG_PATH as needed.
CONFIG_PATH = os.path.join(os.path.dirname(__file__), "config.yaml")
with open(CONFIG_PATH, 'r') as f:
    config = yaml.safe_load(f)

paths_config = config['paths']
EvoEF2_dir = paths_config['EvoEF2_dir']


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
    return parser.parse_args()

def safe_makedir(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)


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
    iface_array = np.array(iface_array)
    while (iface_array.max() - iface_array.min()) > (length - 1):
        if iface_array.std() == 0:
            break
        zscores = np.abs((iface_array - iface_array.mean()) / iface_array.std())
        idx = np.argmax(zscores)
        iface_array = np.delete(iface_array, idx)
    return int(iface_array.min()), int(iface_array.max())

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
    cwd = os.path.abspath(EvoEF2_dir)
    EvoEF2_exec = os.path.join(os.path.abspath(EvoEF2_dir), "EvoEF2")
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


def run_fix_broken_pdb_chain(input_pdb, output_pdb):
    fixed_lines = []
    with open(input_pdb, 'r') as infile:
        for line in infile:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and len(line) >= 78:
                record     = line[0:6]
                serial     = line[6:11]
                spacer1    = line[11:12]
                atom_name  = line[12:16]
                altLoc     = line[16:17]
                resName    = line[17:20]
                spacer2    = line[20:21]
                chain      = line[21:22]
                resSeq     = line[22:26]
                iCode      = line[26:27]
                spacer3    = line[27:30]
                x          = line[30:38]
                y          = line[38:46]
                z          = line[46:54]
                occupancy  = line[54:60]
                tempFactor = line[60:66]
                authChain_field = line[66:76]
                element    = line[76:78]
                remainder  = line[78:]

                new_authChain = f"{chain:>{len(authChain_field)}}"

                new_line = (
                    record + serial + spacer1 + atom_name + altLoc +
                    resName + spacer2 + chain + resSeq + iCode + spacer3 +
                    x + y + z + occupancy + tempFactor +
                    new_authChain + element + remainder
                )
                fixed_lines.append(new_line)
            else:
                fixed_lines.append(line)

    with open(output_pdb, 'w') as outfile:
        outfile.writelines(fixed_lines)


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
        # If out_dir is provided, copy the fixed pdb file to out_dir with new filename.
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
        # Clean up if requested.
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

    iface_array = parse_iface(args.iface)
    min_iface, max_iface = remove_outliers(iface_array, args.length)

    trunc_resi_start = list(range(max_iface - args.length + 1, min_iface + 1))
    trunc_resi_end = [start + args.length - 1 for start in trunc_resi_start]

    summary = []
    for start, end in zip(trunc_resi_start, trunc_resi_end):
        out_filename = f"{input_basename}_{start}_{end}.pdb"
        out_pdb = os.path.join(trunc_dir, out_filename)
        truncate_pdb_by_residue(fixed_pdb, out_pdb, start, end)
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
            "EvoEF2": evoef2_total
        })
        if args.verbose:
            print(f"Processed {out_pdb} with EvoEF2 score {evoef2_total}")

    csv_path = os.path.join(base_temp_dir, f"{input_basename}_EvoEF2_trunc.csv")
    with open(csv_path, "w", newline="") as csvfile:
        fieldnames = ["pdb_path", "trunc_resi_start", "trunc_resi_end", "EvoEF2"]
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

    # If out_dir is provided, copy the best pdb file to that directory with the updated filename.
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

    # If --clean flag is set, remove the temporary directory with intermediate files.
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
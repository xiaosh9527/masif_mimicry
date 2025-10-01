#!/usr/bin/env python3
import argparse
import sys

def fix_pdb(input_path, output_path):
    fixed_lines = []
    with open(input_path, 'r') as infile:
        for line in infile:
            # Only process ATOM/HETATM records that are long enough
            if (line.startswith("ATOM") or line.startswith("HETATM")) and len(line) >= 78:
                # PDB columns (using 0-index slicing):
                #  0-6:   record name
                #  6-11:  atom serial
                # 11-12:  blank
                # 12-16:  atom name
                # 16-17:  altLoc
                # 17-20:  resName
                # 20-21:  blank
                # 21-22:  chain id
                # 22-26:  resSeq
                # 26-27:  iCode
                # 27-30:  blank
                # 30-38:  x
                # 38-46:  y
                # 46-54:  z
                # 54-60:  occupancy
                # 60-66:  tempFactor
                # 66-76:  auth chain id (bad, to be replaced)
                # 76-78:  element
                # 78-end: charge and/or trailing text
                record    = line[0:6]
                serial    = line[6:11]
                spacer1   = line[11:12]
                atom_name = line[12:16]
                altLoc    = line[16:17]
                resName   = line[17:20]
                spacer2   = line[20:21]
                chain     = line[21:22]
                resSeq    = line[22:26]
                iCode     = line[26:27]
                spacer3   = line[27:30]
                x         = line[30:38]
                y         = line[38:46]
                z         = line[46:54]
                occupancy = line[54:60]
                tempFactor= line[60:66]
                # The auth chain id field occupies columns 66-76.
                authChain_field = line[66:76]
                element   = line[76:78]
                remainder = line[78:]
                
                # Replace the auth chain id with the chain id (right-justified in same width)
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
                
    with open(output_path, 'w') as outfile:
        outfile.writelines(fixed_lines)

def main():
    parser = argparse.ArgumentParser(
        description="Fix PDB files by replacing the auth chain id field with the correct chain id."
    )
    parser.add_argument("--input", required=True, help="Path to the input pdb file")
    parser.add_argument("--output", help="Path to the output pdb file; if not given, the input file is overwritten")
    args = parser.parse_args()

    input_path = args.input
    output_path = args.output if args.output else args.input

    try:
        fix_pdb(input_path, output_path)
    except Exception as e:
        sys.exit(f"Error processing file: {e}")

if __name__ == '__main__':
    main()
"""CLI entry point: python -m masif_mimicry.postprocess"""

import os
import sys
from argparse import ArgumentParser
from pathlib import Path

from masif_mimicry.config.paths import resolve_database_paths
from masif_mimicry.postprocess.discovery import discover_deduplicated_rows
from masif_mimicry.postprocess.pipeline import process_results_mimicry
from datetime import datetime


def main(argv=None):
    parser = ArgumentParser(description="Post-process mimicry search CSV outputs")
    parser.add_argument(
        "--target_run_dir",
        type=Path,
        required=True,
        help="Search run directory containing per-P1 subfolders and CSVs",
    )
    parser.add_argument(
        "--target_pdb",
        type=Path,
        required=True,
        help="Path to target chain-C PDB (e.g. 021structure_C.pdb)",
    )
    parser.add_argument(
        "--database_dir",
        type=Path,
        required=True,
        help="MaSIF database root (e.g. TED_domainome/output)",
    )
    parser.add_argument(
        "--target_preprocess_dir",
        type=Path,
        required=True,
        help="MaSIF preprocess directory for the target (e.g. data/NUP98)",
    )
    parser.add_argument("-o", "--out_csv_file", type=Path, required=True, help="Combined output CSV")
    parser.add_argument("--subset", type=Path, default=None, help="Optional P1_id list (one per line)")
    parser.add_argument("--ligand", type=str, required=True, help="Ligand CHAIN_RESNAME, e.g. 'A_021'")

    args = parser.parse_args(argv)

    database_dir = os.path.abspath(os.path.expanduser(str(args.database_dir)))
    _, db_prep = resolve_database_paths(database_dir)
    if not os.path.isdir(db_prep):
        print(f"Error: data_preparation directory not found: {db_prep}", file=sys.stderr)
        sys.exit(1)

    df = discover_deduplicated_rows(args.target_run_dir)
    if df.empty:
        print("No hits to post-process.", flush=True)
        sys.exit(0)

    if args.subset is not None:
        with open(args.subset) as f:
            subset_list = {line.strip() for line in f if line.strip()}
        df = df[df["P1_id"].isin(subset_list)]

    print(f"Starting post-processing at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    process_results_mimicry(
        df,
        target_pdb=args.target_pdb,
        database_dir=database_dir,
        target_preprocess_dir=args.target_preprocess_dir,
        ligand_def=args.ligand,
        out_csv_file=args.out_csv_file,
    )
    print(f"Results written to {args.out_csv_file}")
    print(f"Post-processing completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()

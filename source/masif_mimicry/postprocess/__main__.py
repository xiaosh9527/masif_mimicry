"""CLI entry point: python -m masif_mimicry.postprocess"""

import os
import sys
from argparse import ArgumentParser
from pathlib import Path

from masif_mimicry.config.paths import resolve_database_paths
from masif_mimicry.postprocess.discovery import (
    discover_deduplicated_rows,
    discover_deduplicated_rows_for_subset,
)
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
    parser.add_argument(
        "--cluster_method",
        choices=["max_score", "medoid"],
        default="max_score",
        help="Representative selection within each cluster during deduplication: "
        "max_score (highest MaSIF-score) or medoid (most geometrically central pose, "
        "lowest cluster_mean_rmsd). Affects the kept pose and all downstream metrics.",
    )
    parser.add_argument("--subset", type=Path, default=None, help="Optional P1_id list (one per line)")
    parser.add_argument(
        "--p2_id",
        type=str,
        default=None,
        help="Target P2_id for search CSV paths (required with --subset), e.g. 5B8I_AB_C",
    )
    parser.add_argument(
        "--seed_ppi_id",
        choices=["p1", "p2"],
        default="p1",
        help="Seed PPI id in search CSV filename (default: p1)",
    )
    parser.add_argument(
        "--target_ppi_id",
        choices=["p1", "p2"],
        default="p1",
        help="Target PPI id in search CSV filename (default: p1)",
    )
    parser.add_argument("--ligand", type=str, required=True, help="Ligand CHAIN_RESNAME, e.g. 'A_021'")
    parser.add_argument(
        "--database_info_csv",
        type=Path,
        required=True,
        help="Domain metadata CSV (e.g. TED_human_domainome_info_intracellular.csv); "
        "only rows whose id is in the postprocess P1_id set are loaded into memory.",
    )

    args = parser.parse_args(argv)

    database_info_csv = Path(os.path.abspath(os.path.expanduser(str(args.database_info_csv))))
    if not database_info_csv.is_file():
        print(f"Error: database_info_csv not found: {database_info_csv}", file=sys.stderr)
        sys.exit(1)

    database_dir = os.path.abspath(os.path.expanduser(str(args.database_dir)))
    _, db_prep = resolve_database_paths(database_dir)
    if not os.path.isdir(db_prep):
        print(f"Error: data_preparation directory not found: {db_prep}", file=sys.stderr)
        sys.exit(1)

    if args.subset is not None:
        if args.p2_id is None:
            print("Error: --p2_id is required when --subset is given", file=sys.stderr)
            sys.exit(1)
        with open(args.subset) as f:
            p1_ids = [line.strip() for line in f if line.strip()]
        df = discover_deduplicated_rows_for_subset(
            args.target_run_dir,
            p1_ids,
            args.p2_id,
            seed_ppi_id=args.seed_ppi_id,
            target_ppi_id=args.target_ppi_id,
            database_dir=database_dir,
            target_preprocess_dir=args.target_preprocess_dir,
            cluster_method=args.cluster_method,
        )
    else:
        df = discover_deduplicated_rows(
            args.target_run_dir,
            database_dir=database_dir,
            target_preprocess_dir=args.target_preprocess_dir,
            cluster_method=args.cluster_method,
        )

    if df.empty:
        print("No hits to post-process.", flush=True)
        sys.exit(0)

    print(f"Starting post-processing at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    process_results_mimicry(
        df,
        target_pdb=args.target_pdb,
        database_dir=database_dir,
        target_preprocess_dir=args.target_preprocess_dir,
        ligand_def=args.ligand,
        out_csv_file=args.out_csv_file,
        database_info_csv=database_info_csv,
    )
    print(f"Results written to {args.out_csv_file}")
    print(f"Post-processing completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()

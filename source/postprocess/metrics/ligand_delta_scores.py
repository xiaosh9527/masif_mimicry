#!/usr/bin/env python
"""
Compute ligand-dependent delta scores for one (target, match) pair.

Workflow:
1) Ensure target-no-ligand preprocessing exists (or run preprocess_pdb.py).
2) Run reverse score-only twice via masif_search.py --score_binder --transform:
   - binder = target
   - binder = target-no-ligand
3) Compute deltas: with_ligand - no_ligand.
"""

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


def parse_flattened_transform(flattened_transform: str) -> np.ndarray:
    parts = [x.strip() for x in flattened_transform.split(",") if x.strip()]
    if len(parts) != 16:
        raise ValueError(
            "flattened_transform must contain exactly 16 comma-separated values."
        )
    try:
        mat = np.asarray([float(x) for x in parts], dtype=float).reshape(4, 4)
    except ValueError as exc:
        raise ValueError("flattened_transform contains non-numeric values.") from exc
    if not np.all(np.isfinite(mat)):
        raise ValueError("flattened_transform contains non-finite values.")
    return mat


def invert_transform(flattened_transform: str) -> str:
    transform = parse_flattened_transform(flattened_transform)
    inv_transform = np.linalg.inv(transform)
    return ",".join(str(x) for x in inv_transform.flatten())


def parse_score_line(stdout: str) -> Dict[str, object]:
    lines = [line for line in stdout.splitlines() if line.startswith("query_name:")]
    if not lines:
        raise ValueError("Could not find score output line starting with 'query_name:'.")
    score_line = lines[-1]
    score_dict: Dict[str, object] = {}
    for item in score_line.split(", "):
        if ": " not in item:
            continue
        key, value = item.split(": ", 1)
        try:
            if any(c in value for c in [".", "e", "E"]):
                score_dict[key] = float(value)
            else:
                score_dict[key] = int(value)
        except ValueError:
            score_dict[key] = value
    if "nn_score" not in score_dict:
        raise ValueError("Parsed score line is missing 'nn_score'.")
    return score_dict


def run_score_binder(
    query_preproc_dir: Path,
    query_name: str,
    query_vix: int,
    binder_preproc_dir: Path,
    binder_name: str,
    flattened_transform: str,
    masif_neosurf_repo: Path,
    verbose: bool = False,
) -> Dict[str, object]:
    cmd = [
        sys.executable,
        str(masif_neosurf_repo / "masif_search.py"),
        "--database",
        str(binder_preproc_dir),
        "--score_binder",
        str(binder_name),
        f"--transform={flattened_transform}",
        # "--transform",
        # str(flattened_transform),
        "--target_dir",
        str(query_preproc_dir),
        "--target",
        str(query_name),
        "--site_vix",
        str(query_vix),
    ]
    if verbose:
        print("[run_score_binder] Command:")
        print(" ".join(cmd))

    res = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    if verbose and res.stdout:
        print("[run_score_binder] STDOUT:")
        print(res.stdout)
    if verbose and res.stderr:
        print("[run_score_binder] STDERR:")
        print(res.stderr, file=sys.stderr)
    if res.returncode != 0:
        raise RuntimeError(
            "masif_search.py failed.\n"
            f"Command: {' '.join(cmd)}\n"
            f"STDOUT:\n{res.stdout}\n"
            f"STDERR:\n{res.stderr}"
        )
    return parse_score_line(res.stdout)


def get_target_no_ligand_name(target_name: str) -> str:
    return f"{target_name}-no-ligand"


def get_target_no_ligand_required_paths(
    target_preproc_dir: Path,
    target_name: str,
) -> List[Path]:
    target_no_ligand = get_target_no_ligand_name(target_name)
    return [
        target_preproc_dir
        / "descriptors/sc05/all_feat"
        / target_no_ligand
        / "p1_desc_flipped.npy",
        target_preproc_dir
        / "descriptors/sc05/all_feat"
        / target_no_ligand
        / "p1_desc_straight.npy",
    ]


def has_target_no_ligand(target_preproc_dir: Path, target_name: str) -> bool:
    return all(path.exists() for path in get_target_no_ligand_required_paths(target_preproc_dir, target_name))


def ensure_target_no_ligand(
    target_preproc_dir: Path,
    target_name: str,
    masif_neosurf_repo: Path,
    verbose: bool = False,
) -> Tuple[str, bool]:
    target_no_ligand = get_target_no_ligand_name(target_name)
    if has_target_no_ligand(target_preproc_dir, target_name):
        return target_no_ligand, False

    target_pdb = (
        target_preproc_dir
        / "data_preparation/01-benchmark_pdbs"
        / f"{target_name}.pdb"
    )
    if not target_pdb.exists():
        raise FileNotFoundError(
            f"Target input PDB not found for preprocessing: {target_pdb}"
        )

    cmd = [
        sys.executable,
        str(masif_neosurf_repo / "preprocess_pdb.py"),
        str(target_pdb),
        target_no_ligand,
        "-o",
        str(target_preproc_dir),
    ]
    if verbose:
        print("[ensure_target_no_ligand] Missing no-ligand descriptors; running preprocess.")
        print("[ensure_target_no_ligand] Command:")
        print(" ".join(cmd))
    res = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    if verbose and res.stdout:
        print("[ensure_target_no_ligand] STDOUT:")
        print(res.stdout)
    if verbose and res.stderr:
        print("[ensure_target_no_ligand] STDERR:")
        print(res.stderr, file=sys.stderr)
    if res.returncode != 0:
        raise RuntimeError(
            "preprocess_pdb.py failed while preparing target-no-ligand.\n"
            f"Command: {' '.join(cmd)}\n"
            f"STDOUT:\n{res.stdout}\n"
            f"STDERR:\n{res.stderr}"
        )
    return target_no_ligand, True


def compute_delta_scores(
    target_name: str,
    match_name: str,
    match_vix: int,
    flattened_transform: str,
    target_preproc_dir: Path,
    match_preproc_dir: Path,
    masif_neosurf_repo: Path,
    target_no_ligand_name: str,
    verbose: bool = False,
) -> Dict[str, object]:
    inv_transform = invert_transform(flattened_transform)

    rev_with_ligand = run_score_binder(
        query_preproc_dir=match_preproc_dir,
        query_name=match_name,
        query_vix=match_vix,
        binder_preproc_dir=target_preproc_dir,
        binder_name=target_name,
        flattened_transform=inv_transform,
        masif_neosurf_repo=masif_neosurf_repo,
        verbose=verbose,
    )
    rev_no_ligand = run_score_binder(
        query_preproc_dir=match_preproc_dir,
        query_name=match_name,
        query_vix=match_vix,
        binder_preproc_dir=target_preproc_dir,
        binder_name=target_no_ligand_name,
        flattened_transform=inv_transform,
        masif_neosurf_repo=masif_neosurf_repo,
        verbose=verbose,
    )

    result: Dict[str, object] = {
        "target": target_name,
        "target_no_ligand": target_no_ligand_name,
        "matched_protein": match_name,
        "matched_vix": int(match_vix),
    }
    for key, value in rev_with_ligand.items():
        result[f"rev_with_ligand_{key}"] = value
    for key, value in rev_no_ligand.items():
        result[f"rev_no_ligand_{key}"] = value

    result["delta_nn_score"] = float(rev_with_ligand["nn_score"]) - float(
        rev_no_ligand["nn_score"]
    )
    result["delta_desc_dist_score"] = float(
        rev_with_ligand["desc_dist_score"]
    ) - float(rev_no_ligand["desc_dist_score"])
    result["delta_mean_desc_dist_score"] = float(
        rev_with_ligand["mean_desc_dist_score"]
    ) - float(rev_no_ligand["mean_desc_dist_score"])
    return result


def write_one_row_csv(output_csv: Path, result: Dict[str, object]) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(result.keys())
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(result)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compute ligand delta scores for one target/match pair using reverse "
            "score-binder runs."
        )
    )
    parser.add_argument("--target_preproc_dir", type=Path, required=True)
    parser.add_argument("--target_name", type=str, required=True)
    parser.add_argument("--match_preproc_dir", type=Path, required=True)
    parser.add_argument("--match_name", type=str, required=True)
    parser.add_argument("--match_vix", type=int, required=True)
    parser.add_argument("--flattened_transform", type=str, required=True)
    parser.add_argument("--masif-neosurf-repo", dest="masif_neosurf_repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=None, help="Optional output path for one-row CSV.")
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print subprocess commands and their stdout/stderr for debugging.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    masif_neosurf_repo = Path(__file__).resolve().parent

    if args.match_vix < 0:
        parser.error("--match_vix must be non-negative.")

    for dirname, path in [
        ("target_preproc_dir", args.target_preproc_dir),
        ("match_preproc_dir", args.match_preproc_dir),
        ("masif_neosurf_repo", args.masif_neosurf_repo),
    ]:
        if not path.exists():
            parser.error(f"{dirname} does not exist: {path}")

    masif_search_path = args.masif_neosurf_repo / "masif_search.py"
    preprocess_path = args.masif_neosurf_repo / "preprocess_pdb.py"
    if not masif_search_path.exists():
        parser.error(f"Could not find masif_search.py at: {masif_search_path}")
    if not preprocess_path.exists():
        parser.error(f"Could not find preprocess_pdb.py at: {preprocess_path}")

    try:
        parse_flattened_transform(args.flattened_transform)
    except ValueError as exc:
        parser.error(str(exc))

    try:
        target_no_ligand, preprocessed = ensure_target_no_ligand(
            target_preproc_dir=args.target_preproc_dir,
            target_name=args.target_name,
            masif_neosurf_repo=args.masif_neosurf_repo,
            verbose=args.verbose,
        )
        delta_result = compute_delta_scores(
            target_name=args.target_name,
            match_name=args.match_name,
            match_vix=args.match_vix,
            flattened_transform=args.flattened_transform,
            target_preproc_dir=args.target_preproc_dir,
            match_preproc_dir=args.match_preproc_dir,
            masif_neosurf_repo=args.masif_neosurf_repo,
            target_no_ligand_name=target_no_ligand,
            verbose=args.verbose,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    delta_result["no_ligand_preprocessed_now"] = preprocessed
    print(json.dumps(delta_result, indent=2, sort_keys=True))

    if args.output is not None:
        try:
            write_one_row_csv(args.output, delta_result)
        except Exception as exc:
            print(f"ERROR writing output CSV: {exc}", file=sys.stderr)
            return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

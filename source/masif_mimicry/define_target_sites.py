"""Define target surface sites and write a shared target run directory for parallel search."""

import argparse
import os
import shutil
import sys

import numpy as np

from masif_mimicry.config.paths import set_params
from masif_mimicry.search.features import get_features
from masif_mimicry.search.manifest import TARGET_SITES_MANIFEST
from masif_mimicry.search.target_sites import (
    select_target_sites_by_radius_fps,
    write_partner_pdb_for_clashes,
    write_target_sites_manifest,
    write_target_vert_files,
)


def create_parser():
    p = argparse.ArgumentParser(
        description="Select target sites and prepare target_vert/ for masif_mimicry search",
    )
    p.add_argument(
        "--target_preprocess_dir",
        type=str,
        required=True,
        help="MaSIF target preprocess directory (same layout as search)",
    )
    p.add_argument("--target_pdb", type=str, required=True, help="Target PDB id (e.g. 021structure_C_AB)")
    p.add_argument("--target_ppi_id", choices=["p1", "p2"], default="p1")
    p.add_argument("--target_chain", type=str, required=True, help="Target chain id for sanity checks")
    p.add_argument("--target_residue", type=int, required=True, help="Target residue number")
    p.add_argument("--target_atom", type=str, required=True, help="Target atom name (e.g. CB)")
    p.add_argument(
        "--target_sampling_radius",
        type=float,
        default=5.0,
        help="Radius (A) around target atom for candidate vertices",
    )
    p.add_argument(
        "--num_points",
        type=int,
        default=5,
        help="Number of target sites after FPS subsampling",
    )
    p.add_argument(
        "--target_run_dir",
        type=str,
        required=True,
        help="Output directory for this target run (e.g. .../search_results/021structure_C_AB_)",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help=f"Replace existing {TARGET_SITES_MANIFEST} and target_vert/",
    )
    return p


def main(args):
    P2 = args.target_pdb
    target_run_dir = os.path.abspath(os.path.expanduser(args.target_run_dir))

    if args.target_ppi_id == "p1":
        assert args.target_chain == P2.split("_")[1], (
            f"Target chain {args.target_chain} does not match PDB {P2}"
        )
    else:
        assert args.target_chain == P2.split("_")[2], (
            f"Target chain {args.target_chain} does not match PDB {P2}"
        )

    manifest_path = os.path.join(target_run_dir, TARGET_SITES_MANIFEST)
    if os.path.exists(manifest_path) and not args.overwrite:
        print(
            f"Error: {manifest_path} already exists. Use --overwrite to replace.",
            file=sys.stderr,
        )
        sys.exit(1)

    target_preprocess_dir = os.path.abspath(os.path.expanduser(args.target_preprocess_dir))
    if not os.path.isdir(target_preprocess_dir):
        print(f"Error: target_preprocess_dir not found: {target_preprocess_dir}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(target_run_dir, exist_ok=True)
    params = set_params(database_dir=target_preprocess_dir, target_preprocess_dir=target_preprocess_dir)
    P2_all_feats = get_features(params, P2, args.target_ppi_id, source=False, flip_desc=False)

    try:
        selected_points_idx = select_target_sites_by_radius_fps(
            np.array(P2_all_feats["mesh"].vertices),
            P2_all_feats["pdb"],
            chain=args.target_chain,
            residue=args.target_residue,
            atom_name=args.target_atom,
            radius=args.target_sampling_radius,
            num_points=args.num_points,
        )
    except ValueError as e:
        print(f"Error: {e}", flush=True)
        sys.exit(1)

    sites = [int(x) for x in np.asarray(selected_points_idx).tolist()]
    target_pdb_filename = f"{P2.split('_')[0]}_{args.target_chain}.pdb"
    partner_pdb_filename, _ = write_partner_pdb_for_clashes(
        params, P2, args.target_ppi_id, target_run_dir,
    )

    write_target_vert_files(P2_all_feats, selected_points_idx, args.target_ppi_id, target_run_dir)
    shutil.copy(
        P2_all_feats["pdb"],
        os.path.join(target_run_dir, target_pdb_filename),
    )

    manifest = {
        "target_pdb": P2,
        "target_ppi_id": args.target_ppi_id,
        "target_chain": args.target_chain,
        "selection": "residue_fps",
        "target_residue": args.target_residue,
        "target_atom": args.target_atom,
        "target_sampling_radius": args.target_sampling_radius,
        "num_points": args.num_points,
        "sites": sites,
        "target_pdb_file": target_pdb_filename,
        "partner_pdb_file": partner_pdb_filename,
    }
    write_target_sites_manifest(target_run_dir, manifest)

    print(
        f"Prepared {len(sites)} target sites in {target_run_dir}\n"
        f"  manifest: {manifest_path}\n"
        f"  target_vert/: {len(sites)} patch files\n"
        f"  {target_pdb_filename}, {partner_pdb_filename}",
        flush=True,
    )


if __name__ == "__main__":
    main(create_parser().parse_args())

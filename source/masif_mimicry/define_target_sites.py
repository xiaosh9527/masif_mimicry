"""Define target surface sites and write a shared target run directory for parallel search."""

import argparse
import os
import shutil
import sys

import numpy as np

from masif_mimicry.config.paths import set_params
from masif_mimicry.search.features import get_features
from masif_mimicry.search.manifest import TARGET_SITES_MANIFEST
from masif_mimicry.search.docking import select_patches
from masif_mimicry.search.target_sites import (
    parse_site_idx_list,
    select_target_sites_by_grid,
    validate_site_indices,
    write_partner_pdb_for_clashes,
    write_target_sites_manifest,
    write_target_vert_files,
)
from masif_mimicry.utils.parse_structures import target_chain_from_pdb_id


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
    p.add_argument(
        "--query_pdb",
        type=str,
        help="PDB file containing the query residue (any path; same frame as preprocessed target)",
    )
    p.add_argument(
        "--query_chain",
        type=str,
        help="Chain id in --query_pdb for the query residue",
    )
    p.add_argument(
        "--query_residue",
        type=int,
        help="Residue number in --query_pdb for grid site selection",
    )
    p.add_argument(
        "--grid_distance_cutoff",
        type=float,
        default=4.0,
        help="Distance (A) from any query-residue heavy atom to include mesh vertices",
    )
    p.add_argument(
        "--num_points",
        type=int,
        default=5,
        help="Number of target sites (grid mode) or minimum when using exhaustive selection",
    )
    p.add_argument("--iface_cutoff", type=float, default=0.0, help="Interface score cutoff (exhaustive mode)")
    p.add_argument(
        "--top_iface_percent",
        type=float,
        default=0.0,
        help="Top iface percentile cap (exhaustive mode)",
    )
    p.add_argument(
        "--interface_only",
        action="store_true",
        help="Exhaustive mode: only interface-labeled vertices",
    )
    p.add_argument(
        "--downsample",
        type=int,
        default=1,
        help="Downsample rate for exhaustive site selection only",
    )
    p.add_argument(
        "--site_idx",
        type=str,
        default=None,
        help="Comma-separated MaSIF patch center indices (e.g. '1207,940,300'). "
        "Overrides grid/exhaustive selection when set.",
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

    expected_chain = target_chain_from_pdb_id(P2, args.target_ppi_id)
    assert args.target_chain == expected_chain, (
        f"Target chain {args.target_chain} does not match PDB {P2} (expected {expected_chain})"
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
    n_patch_centers = len(P2_all_feats["desc"])

    if args.site_idx is not None:
        grid_args = (args.query_pdb, args.query_chain, args.query_residue)
        if any(x is not None for x in grid_args):
            print(
                "WARNING: --site_idx overrides --query_pdb/--query_chain/--query_residue.",
                flush=True,
            )
        try:
            site_list = parse_site_idx_list(args.site_idx)
            selected_points_idx = validate_site_indices(site_list, n_patch_centers, P2)
        except ValueError as e:
            print(f"Error: {e}", flush=True)
            sys.exit(1)
        selection = "manual"
        print(
            f"Using {len(selected_points_idx)} manual site indices: {selected_points_idx.tolist()}",
            flush=True,
        )
    elif args.query_pdb and args.query_chain and args.query_residue is not None:
        query_pdb = os.path.abspath(os.path.expanduser(args.query_pdb))
        if not os.path.isfile(query_pdb):
            print(f"Error: query_pdb not found: {query_pdb}", file=sys.stderr)
            sys.exit(1)
        try:
            selected_points_idx = select_target_sites_by_grid(
                P2_all_feats,
                query_pdb=query_pdb,
                chain=args.query_chain,
                residue=args.query_residue,
                num_points=args.num_points,
                distance_cutoff=args.grid_distance_cutoff,
            )
        except ValueError as e:
            print(f"Error: {e}", flush=True)
            sys.exit(1)
        selection = "grid"
        print(
            f"Selected {len(selected_points_idx)} grid sites near residue "
            f"{args.query_chain}:{args.query_residue} from {query_pdb} "
            f"(cutoff={args.grid_distance_cutoff} A, num_points={args.num_points}).",
            flush=True,
        )
    elif args.query_pdb or args.query_chain or args.query_residue is not None:
        print(
            "Error: grid selection requires --query_pdb, --query_chain, and --query_residue together.",
            file=sys.stderr,
        )
        sys.exit(1)
    else:
        selected_points_idx, _, _ = select_patches(
            P2_all_feats,
            downsample_rate=args.downsample,
            iface_cutoff=args.iface_cutoff,
            top_iface_percent=args.top_iface_percent,
            interface_only=args.interface_only,
            verbose=True,
        )
        if len(selected_points_idx) == 0:
            print("Error: No target sites selected. Check --downsample and iface flags.", file=sys.stderr)
            sys.exit(1)
        selection = "exhaustive"

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
        "selection": selection,
        "num_points": args.num_points,
        "sites": sites,
        "target_pdb_file": target_pdb_filename,
        "partner_pdb_file": partner_pdb_filename,
    }
    if selection == "manual":
        manifest["site_idx"] = args.site_idx
    if selection == "grid":
        manifest["query_pdb"] = os.path.abspath(os.path.expanduser(args.query_pdb))
        manifest["query_chain"] = args.query_chain
        manifest["query_residue"] = args.query_residue
        manifest["grid_distance_cutoff"] = args.grid_distance_cutoff
    if selection == "exhaustive":
        manifest["downsample"] = args.downsample
        manifest["iface_cutoff"] = args.iface_cutoff
        manifest["top_iface_percent"] = args.top_iface_percent
        manifest["interface_only"] = args.interface_only
    write_target_sites_manifest(target_run_dir, manifest)

    partner_note = partner_pdb_filename if partner_pdb_filename else "(no partner chains)"
    print(
        f"Prepared {len(sites)} target sites in {target_run_dir}\n"
        f"  manifest: {manifest_path}\n"
        f"  target_vert/: {len(sites)} patch files\n"
        f"  {target_pdb_filename}, partner: {partner_note}",
        flush=True,
    )


if __name__ == "__main__":
    main(create_parser().parse_args())

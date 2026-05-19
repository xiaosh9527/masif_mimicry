import os, sys, argparse, shutil, tempfile
import numpy as np
import pandas as pd

from utils import *
from subprocess import Popen, PIPE

# Target-site sampling flags (use define_target_sites.py when --target_run_dir is set).
TARGET_SITE_ARG_NAMES = (
    "target_residue",
    "target_atom",
    "target_sampling_radius",
    "num_points",
)


def create_parser():
    p = argparse.ArgumentParser("Simplified MaSIF mimicry search")

    p.add_argument("--database_dir", type=str, required=True,
                   help="MaSIF seed database directory (preprocessed data/<db_name>/ tree)")
    p.add_argument("--target_preprocess_dir", type=str, required=True,
                   help="MaSIF target preprocess directory (same layout for the target protein)")

    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--seed_pdb", type=str, help="Single seed PDB (format: PDB_ppi_chain or PDB_ppi)")
    group.add_argument("--split_seed_list", type=str, help="File containing one seed PDB id per line")

    p.add_argument(
        "--target_run_dir",
        type=str,
        default=None,
        help="Prepared target run directory (target_sites.json + target_vert/). "
             "When set, target site selection flags are not allowed.",
    )

    p.add_argument("--target_pdb", type=str, default=None,
                   help="Target PDB identifier (legacy mode without --target_run_dir)")
    p.add_argument("--target_ppi_id", choices=["p1", "p2"], default="p1",
                   help="ppi side of the target (legacy mode only; read from manifest if --target_run_dir)")
    p.add_argument("--target_chain", type=str, default=None,
                   help="Target chain id (legacy mode only; read from manifest if --target_run_dir)")

    p.add_argument("--target_residue", type=int, default=None,
                   help="[define_target_sites.py only] Target residue number")
    p.add_argument("--target_atom", type=str, default=None,
                   help="[define_target_sites.py only] Target atom name")
    p.add_argument("--target_sampling_radius", type=float, default=None,
                   help="[define_target_sites.py only] Radius (A) around target atom")
    p.add_argument("--num_points", type=int, default=None,
                   help="[define_target_sites.py only] FPS subsample count")

    p.add_argument("--downsample", type=int, default=1, help="Downsample rate for seed site selection")
    p.add_argument("--top_iface_percent", type=float, default=0.0,
                   help="Top percentile of seed patches by interface score")
    p.add_argument("--iface_cutoff", type=float, default=0.0, help="Interface score cutoff for seed patches")
    p.add_argument("--interface_only", action="store_true",
                   help="Search only seed interfaces (requires seed complexes)")
    p.add_argument("--desc_dist_cutoff", type=float, default=1.5, help="Descriptor distance cutoff (filtering)")
    p.add_argument("--desc_dist_score_cutoff", type=float, default=0.25, help="Post-alignment score cutoff")

    p.add_argument("--output_dir", type=str, default=".", help="Output directory (legacy mode only)")
    p.add_argument("--output_postfix", type=str, default="", help="Subfolder postfix (legacy mode only)")
    p.add_argument("--count_clashes", action="store_true", help="Compute clashes (slower)")
    p.add_argument("--ca_clash_threshold", type=float, default=100, help="CA clash threshold")
    p.add_argument("--heavy_atom_clash_threshold", type=float, default=100, help="Heavy-atom clash threshold")
    p.add_argument("--compute_source_residues", action="store_true",
                   help="Map patch centers to residues via surf2atom (DSSP; slower)")

    return p


def log(msg):
    print(msg, flush=True)


def flatten_transform(T, precision=16):
    """Serialize a 4x4 transform as comma-separated floats (row-major)."""
    return ",".join(f"{x:.{precision}f}" for x in np.asarray(T).reshape(-1))


def get_usalign_tmscores(p1_pdb, p2_pdb, cache):
    """Run USalign once per (p1_pdb, p2_pdb) pair; return (TMscore_P1, TMscore_P2)."""
    key = (os.path.abspath(p1_pdb), os.path.abspath(p2_pdb))
    if key not in cache:
        process = Popen(
            ["/install/USalign/USalign", p1_pdb, p2_pdb, "-mm", "0", "-ter", "2"],
            stdout=PIPE, stderr=PIPE,
        )
        stdout, _ = process.communicate()
        tmscore_list = [
            float(x.split(" ")[1])
            for x in stdout.decode().splitlines()
            if "TM-score=" in x
        ]
        cache[key] = (
            tmscore_list[0] if tmscore_list else 0.0,
            tmscore_list[1] if len(tmscore_list) > 1 else 0.0,
        )
    return cache[key]


def _reject_target_site_args(args):
    if not args.target_run_dir:
        return
    bad = [f"--{name}" for name in TARGET_SITE_ARG_NAMES if getattr(args, name) is not None]
    if bad:
        print(
            "Error: target site selection flags are not allowed with --target_run_dir:\n  "
            + ", ".join(bad)
            + "\nRun define_target_sites.py to prepare the target run directory.",
            file=sys.stderr,
        )
        sys.exit(1)


def _validate_target_chain(target_ppi_id, target_chain, P2):
    if target_ppi_id == "p1":
        assert target_chain == P2.split("_")[1], f"Target chain {target_chain} does not match PDB {P2}"
    else:
        assert target_chain == P2.split("_")[2], f"Target chain {target_chain} does not match PDB {P2}"


def _resolve_target_setup(args, params):
    """Return (P2, target_ppi_id, target_chain, P2_selected_points_idx, P2_patch_descs, p2_output_root, manifest)."""
    if args.target_run_dir:
        p2_output_root = os.path.abspath(os.path.expanduser(args.target_run_dir))
        manifest, P2_selected_points_idx = load_target_run_manifest(p2_output_root)
        P2 = manifest["target_pdb"]
        target_ppi_id = manifest["target_ppi_id"]
        target_chain = manifest["target_chain"]
        _validate_target_chain(target_ppi_id, target_chain, P2)
        P2_all_feats = get_features(params, P2, target_ppi_id, source=False, flip_desc=False)
        P2_patch_descs = P2_all_feats["desc"][P2_selected_points_idx]
        log(
            f"Using {len(P2_selected_points_idx)} pre-defined target sites from {p2_output_root} "
            f"(residue {manifest['target_residue']} {manifest['target_atom']} chain {target_chain})"
        )
        return P2, target_ppi_id, target_chain, P2_selected_points_idx, P2_patch_descs, p2_output_root, manifest

    if not args.target_pdb or not args.target_chain:
        print("Error: --target_pdb and --target_chain are required without --target_run_dir", file=sys.stderr)
        sys.exit(1)

    P2 = args.target_pdb
    target_ppi_id = args.target_ppi_id
    target_chain = args.target_chain
    _validate_target_chain(target_ppi_id, target_chain, P2)

    P2_all_feats = get_features(params, P2, target_ppi_id, source=False, flip_desc=False)
    radius = args.target_sampling_radius if args.target_sampling_radius is not None else 5.0
    num_points = args.num_points if args.num_points is not None else 5

    if args.target_residue and args.target_atom:
        try:
            P2_selected_points_idx = select_target_sites_by_radius_fps(
                np.array(P2_all_feats["mesh"].vertices),
                P2_all_feats["pdb"],
                chain=target_chain,
                residue=args.target_residue,
                atom_name=args.target_atom,
                radius=radius,
                num_points=num_points,
            )
        except ValueError as e:
            print(f"Error: {e}", flush=True)
            sys.exit(1)
        P2_patch_descs = P2_all_feats["desc"][P2_selected_points_idx]
    else:
        P2_selected_points_idx, P2_patch_descs, _ = select_patches(
            P2_all_feats,
            downsample_rate=args.downsample,
            iface_cutoff=args.iface_cutoff,
            top_iface_percent=args.top_iface_percent,
            interface_only=args.interface_only,
            verbose=False,
        )
        if len(P2_selected_points_idx) == 0:
            print(f"No points selected for {P2}. Please check the parameters. Exiting...")
            sys.exit(1)

    p2_output_root = os.path.join(args.output_dir, f"{P2}_{args.output_postfix}")
    os.makedirs(p2_output_root, exist_ok=True)
    write_target_vert_files(P2_all_feats, P2_selected_points_idx, target_ppi_id, p2_output_root)
    shutil.copy(
        P2_all_feats["pdb"],
        os.path.join(p2_output_root, f"{P2.split('_')[0]}_{target_chain}.pdb"),
    )
    return P2, target_ppi_id, target_chain, P2_selected_points_idx, P2_patch_descs, p2_output_root, None


def main(args):
    _reject_target_site_args(args)

    local_tmp_dir = os.getenv("TMPDIR")
    if not local_tmp_dir or not os.path.isdir(local_tmp_dir):
        local_tmp_dir = tempfile.mkdtemp(prefix="masif_mimicry_")
    else:
        os.makedirs(local_tmp_dir, exist_ok=True)

    if args.split_seed_list:
        lines = open(args.split_seed_list, "r").readlines()
        lines = [line.strip() for line in lines]
        print(f"Number of proteins to search from: {len(lines)}")
    elif args.seed_pdb:
        lines = [args.seed_pdb]
    else:
        print("Please specify either --split_seed_list or --seed_pdb")
        sys.exit(1)

    database_dir = os.path.abspath(os.path.expanduser(args.database_dir))
    target_preprocess_dir = os.path.abspath(os.path.expanduser(args.target_preprocess_dir))
    for label, path in [("database_dir", database_dir), ("target_preprocess_dir", target_preprocess_dir)]:
        if not os.path.isdir(path):
            print(f"Error: {label} does not exist or is not a directory: {path}")
            sys.exit(1)

    params = set_params(database_dir=database_dir, target_preprocess_dir=target_preprocess_dir)
    (
        P2,
        target_ppi_id,
        target_chain,
        P2_selected_points_idx,
        P2_patch_descs,
        p2_output_root,
        manifest,
    ) = _resolve_target_setup(args, params)

    P2_all_feats = get_features(params, P2, target_ppi_id, source=False, flip_desc=False)

    if not args.target_run_dir:
        if args.target_residue and args.target_atom:
            print(
                f"Searching sites similar to {P2} from {','.join(set(lines))}.\n"
                f"This will go through {len(P2_selected_points_idx)} points...",
                flush=True,
            )
        else:
            print(
                f"Searching sites similar to {P2} from {','.join(set(lines))}.\n"
                f"WARNING: Exhaustive alignment mode ({len(P2_selected_points_idx)} points)...",
                flush=True,
            )

    scores = {}
    for P1 in lines:
        scores[(P1, P2)] = {
            "P1_id": [],
            "P2_id": [],
            "P1_source_ppi_id": [],
            "P2_source_ppi_id": [],
            "P1_source_site": [],
            "P2_source_site": [],
            "P1_source_TMscore": [],
            "P2_source_TMscore": [],
            "P1_source_iface": [],
            "P2_source_iface": [],
            "MaSIF-score": [],
            "flattened_transform": [],
        }
        if args.compute_source_residues:
            scores[(P1, P2)]["P1_source_residue"] = []
            scores[(P1, P2)]["P2_source_residue"] = []
        if args.count_clashes:
            scores[(P1, P2)]["ca_clash"] = []
            scores[(P1, P2)]["heavy_atom_clash"] = []

    if not args.interface_only:
        print(f"Aligning to all points with iface > {args.iface_cutoff} and desc dist < {args.desc_dist_cutoff}.")
    else:
        print(f"Aligning to interface points only with desc dist < {args.desc_dist_cutoff}.")

    tmscore_cache = {}
    filtered_target_structure = None
    partner_pdb_path = None
    if args.count_clashes:
        if manifest is not None:
            partner_pdb_path = os.path.join(p2_output_root, manifest["partner_pdb_file"])
            if not os.path.isfile(partner_pdb_path):
                print(f"Error: missing partner PDB {partner_pdb_path}. Re-run define_target_sites.py.", file=sys.stderr)
                sys.exit(1)
            filtered_target_structure = PDBParser(QUIET=True).get_structure("", partner_pdb_path)
        else:
            P2_raw_pdb = os.path.join(
                params["masif_target_root"],
                "data_preparation", "00-raw_pdbs", f"{P2.split('_')[0]}.pdb",
            )
            P2_partner_chain_ids = parse_partner_chain_ids(P2, target_ppi_id)
            P2_partner_chain_suffix = "".join(P2_partner_chain_ids)
            filtered_target_structure = get_filtered_target_structure(
                P2_raw_pdb, P2_partner_chain_ids, {},
            )
            partner_pdb_path = os.path.join(
                p2_output_root, f"{P2.split('_')[0]}_{P2_partner_chain_suffix}.pdb",
            )
            if not os.path.isfile(partner_pdb_path):
                io = PDBIO()
                io.set_structure(filtered_target_structure)
                io.save(partner_pdb_path)

    for P1 in lines:
        if len(P1.split("_")) == 2:
            ppi_id_list = ["p1"]
        elif len(P1.split("_")) == 3:
            ppi_id_list = ["p1", "p2"]
        else:
            raise ValueError("Invalid format. Please use PDB_X or PDB_X_X")

        for ppi_id in ppi_id_list:
            P1_all_feats = get_features(params, P1, ppi_id, source=True, flip_desc=False)
            P1_selected_points_idx, P1_patch_descs, _ = select_patches(
                P1_all_feats,
                downsample_rate=args.downsample,
                iface_cutoff=args.iface_cutoff,
                top_iface_percent=args.top_iface_percent,
                interface_only=args.interface_only,
                verbose=False,
            )

            desc_dist = np.linalg.norm(P1_patch_descs[:, None, :] - P2_patch_descs[None, :, :], axis=-1)

            if np.min(desc_dist) > args.desc_dist_cutoff:
                log(f"Skipping {P1} ({ppi_id}): no descriptors within {args.desc_dist_cutoff}")
                continue

            seed_output_dir = os.path.join(p2_output_root, P1)
            log(f"Aligning {P1} ({ppi_id}) to {P2}...")
            total_matches = 0
            sites_aligned = 0
            need_tmp_pdb = args.count_clashes or args.compute_source_residues

            for i, P2_center in enumerate(P2_selected_points_idx):
                P1_selected_points_idx_final = P1_selected_points_idx[
                    np.where(desc_dist[:, i] < args.desc_dist_cutoff)
                ]

                if len(P1_selected_points_idx_final) == 0:
                    log(f"Skipped site {P2_center} (no descriptors within {args.desc_dist_cutoff})")
                    continue

                sites_aligned += 1
                n_hits = 0

                all_results, _, _, _ = multidock(
                    source_pt=P1_selected_points_idx_final,
                    source_pcd=P1_all_feats["pcd"],
                    source_patch_idxs=P1_all_feats["indices"],
                    source_descs=P1_all_feats["desc"],
                    target_pt=P2_center,
                    target_pcd=P2_all_feats["pcd"],
                    target_patch_idxs=P2_all_feats["indices"],
                    target_descs=P2_all_feats["desc"],
                    binder_align=False,
                )

                for result, P1_center in zip(all_results, P1_selected_points_idx_final):
                    out_filename_base = (
                        f"{P1}_{ppi_id}_{P1_center}_to_{P2}_{target_ppi_id}_{P2_center}"
                    )

                    P1_patch_coords = transform_patch_coords(
                        P1_all_feats["pcd"],
                        P1_all_feats["indices"],
                        P1_center,
                        result.transformation,
                    )
                    descriptor_score = compute_descriptor_score(
                        P1_patch_coords,
                        P1_all_feats["desc"],
                        P1_all_feats["indices"],
                        P1_center,
                        P2_all_feats["pcd"],
                        P2_all_feats["desc"],
                        P2_all_feats["indices"],
                        P2_center,
                    )

                    if descriptor_score < args.desc_dist_score_cutoff:
                        continue

                    tmp_pdb_path = None
                    if need_tmp_pdb:
                        tmp_pdb_path = os.path.join(local_tmp_dir, f"{out_filename_base}.pdb")
                        _ = transform_structure(
                            P1_all_feats["pdb"], result.transformation, tmp_pdb_path,
                        )

                    if args.count_clashes:
                        ca_clashes, heavy_clashes, masif_score = compute_hit_clash_score(
                            tmp_pdb_path,
                            filtered_target_structure,
                            descriptor_score,
                            ca_clash_threshold=args.ca_clash_threshold,
                            heavy_atom_clash_threshold=args.heavy_atom_clash_threshold,
                        )
                        if masif_score < args.desc_dist_score_cutoff:
                            if tmp_pdb_path and os.path.exists(tmp_pdb_path):
                                os.remove(tmp_pdb_path)
                            continue
                    else:
                        ca_clashes, heavy_clashes = 0, 0
                        masif_score = descriptor_score

                    n_hits += 1
                    if n_hits == 1:
                        os.makedirs(seed_output_dir, exist_ok=True)

                    if args.compute_source_residues:
                        target_atom, _ = surf2atom(
                            point_coords=np.array(P2_all_feats["pcd"].points)[P2_center].reshape(1, -1),
                            pdb_path=P2_all_feats["pdb"],
                        )
                        target_residue = target_atom[0].get_parent().get_id()[1]
                        P1_center_coord = apply_transform(
                            np.asarray(P1_all_feats["pcd"].points)[P1_center],
                            result.transformation,
                        )
                        P1_nearest_atom, _ = surf2atom(
                            point_coords=P1_center_coord.reshape(1, -1),
                            pdb_path=tmp_pdb_path,
                        )
                        P1_nearest_res = P1_nearest_atom[0].get_parent().get_id()[1]

                    if tmp_pdb_path and os.path.exists(tmp_pdb_path):
                        os.remove(tmp_pdb_path)

                    TMscore_P1, TMscore_P2 = get_usalign_tmscores(
                        P1_all_feats["pdb"], P2_all_feats["pdb"], tmscore_cache,
                    )

                    scores[(P1, P2)]["P1_id"].append(P1)
                    scores[(P1, P2)]["P2_id"].append(P2)
                    scores[(P1, P2)]["P1_source_ppi_id"].append(ppi_id)
                    scores[(P1, P2)]["P2_source_ppi_id"].append(target_ppi_id)
                    if args.compute_source_residues:
                        scores[(P1, P2)]["P1_source_residue"].append(P1_nearest_res)
                        scores[(P1, P2)]["P2_source_residue"].append(target_residue)
                    scores[(P1, P2)]["P1_source_site"].append(P1_center)
                    scores[(P1, P2)]["P2_source_site"].append(P2_center)
                    scores[(P1, P2)]["P1_source_TMscore"].append(TMscore_P1)
                    scores[(P1, P2)]["P2_source_TMscore"].append(TMscore_P2)
                    scores[(P1, P2)]["P1_source_iface"].append(P1_all_feats["iface"][0][P1_center])
                    scores[(P1, P2)]["P2_source_iface"].append(P2_all_feats["iface"][0][P2_center])
                    scores[(P1, P2)]["MaSIF-score"].append(masif_score)
                    scores[(P1, P2)]["flattened_transform"].append(flatten_transform(result.transformation))
                    if args.count_clashes:
                        scores[(P1, P2)]["ca_clash"].append(ca_clashes)
                        scores[(P1, P2)]["heavy_atom_clash"].append(heavy_clashes)

                total_matches += n_hits
                log(f"Found {n_hits} matches to site {P2_center}")

            log(f"Done {P1} ({ppi_id}): {total_matches} matches across {sites_aligned} sites")

            if len(scores[(P1, P2)]["P1_id"]) > 0:
                hit_df = pd.DataFrame(scores[(P1, P2)])
                hit_df = structural_clusters(hit_df, database_root=params["top_seed_dir"])
                hit_df.to_csv(
                    f"{seed_output_dir}/{P1}_{ppi_id}_to_{P2}_{target_ppi_id}.csv",
                    index=False,
                )


if __name__ == "__main__":
    parser = create_parser()
    main(parser.parse_args())

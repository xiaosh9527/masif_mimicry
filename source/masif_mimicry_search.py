import os, sys, glob, shutil, pickle, csv, argparse, tempfile
import numpy as np
import pandas as pd

from utils import *
from subprocess import Popen, PIPE

def create_parser():
    p = argparse.ArgumentParser("Simplified MaSIF mimicry search")

    # Database / targets
    p.add_argument("--database_dir", type=str, required=True,
                   help="MaSIF seed database directory (preprocessed data/<db_name>/ tree)")
    p.add_argument("--target_preprocess_dir", type=str, required=True,
                   help="MaSIF target preprocess directory (same layout for the target protein)")

    # Which PDBs/chains to use
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--seed_pdb", type=str, help="Single seed PDB (format: PDB_ppi_chain or PDB_ppi)")
    group.add_argument("--split_seed_list", type=str, help="File containing one seed PDB id per line")

    p.add_argument("--target_pdb", type=str, required=True, help="Target PDB identifier (format: PDB_ppi_chain)")
    p.add_argument("--target_ppi_id", choices=["p1", "p2"], default="p1", help="ppi side of the target to align to")
    p.add_argument("--target_chain", type=str, required=True, help="Target chain id (single char) for sanity checks")

    # Optional: allow user to specify a target residue/atom/chain to select surface points near that atom
    p.add_argument("--target_residue", type=int, help="Target residue number (integer). If provided with --target_atom and --target_chain, surface sites are selected within --target_sampling_radius and subsampled with FPS.")
    p.add_argument("--target_atom", type=str, help="Target atom name (e.g., CA, NZ). Used together with --target_residue to find nearby surface vertices.")
    p.add_argument("--target_sampling_radius", type=float, default=5.0,
                   help="Radius (A) around --target_atom for candidate vertices (residue-target mode only)")

    # Selection and thresholds
    p.add_argument("--num_points", type=int, default=5, help="Number of target sites after FPS subsampling (residue-target mode only)")
    p.add_argument("--downsample", type=int, default=1, help="Downsample rate for site selection")
    p.add_argument("--top_iface_percent", type=float, default=0.0, help="Top percentile of points based on interface score value to prioritize")
    p.add_argument("--iface_cutoff", type=float, default=0.0, help="Interface score cutoff")
    p.add_argument("--interface_only", action="store_true", help="Search only seed interfaces (requires seed complexes)")
    p.add_argument("--desc_dist_cutoff", type=float, default=1.5, help="Descriptor distance cutoff (filtering)")
    p.add_argument("--desc_dist_score_cutoff", type=float, default=0.25, help="Post-alignment score cutoff")

    # Output and runtime
    p.add_argument("--output_dir", type=str, default='.', help="Output directory")
    p.add_argument("--output_postfix", type=str, default='', help="Subfolder postfix for outputs")
    p.add_argument("--count_clashes", action='store_true', help="Compute clashes (slower)")
    p.add_argument("--ca_clash_threshold", type=float, default=100, help="CA clash threshold")
    p.add_argument("--heavy_atom_clash_threshold", type=float, default=100, help="Heavy-atom clash threshold")
    p.add_argument("--compute_source_residues", action='store_true',
                   help="Map patch centers to nearest surface-exposed residues via surf2atom (DSSP; slower). "
                        "Adds P1_source_residue and P2_source_residue columns to output CSV.")

    return p


def write_target_vert_files(p2_all_feats, selected_points_idx, target_ppi_id, output_dir, outward_shift=0.25):
    """Write full geodesic patch vertices for each selected target site."""
    vert_dir = os.path.join(output_dir, 'target_vert')
    os.makedirs(vert_dir, exist_ok=True)
    pcd_points = np.asarray(p2_all_feats['pcd'].points)
    centers_path = os.path.join(vert_dir, 'target.vert')
    with open(centers_path, 'w') as out_centers:
        for site_vix in selected_points_idx:
            center = pcd_points[site_vix]
            out_centers.write('{}, {}, {}\n'.format(center[0], center[1], center[2]))
            target_patch, _, _ = get_patch_geo(
                p2_all_feats['pcd'],
                p2_all_feats['indices'],
                site_vix,
                p2_all_feats['desc'],
                flip_normals=False,
                outward_shift=outward_shift,
            )
            vert_path = os.path.join(vert_dir, f'{target_ppi_id}_{site_vix}.vert')
            with open(vert_path, 'w') as out_patch:
                for point in target_patch.points:
                    out_patch.write('{}, {}, {}\n'.format(point[0], point[1], point[2]))
    print(
        f'Wrote {len(selected_points_idx)} patch files and {centers_path} to {vert_dir}',
        flush=True,
    )


def log(msg):
    print(msg, flush=True)


def flatten_transform(T, precision=16):
    """Serialize a 4x4 transform as comma-separated floats (row-major)."""
    return ','.join(f'{x:.{precision}f}' for x in np.asarray(T).reshape(-1))


def get_usalign_tmscores(p1_pdb, p2_pdb, cache):
    """Run USalign once per (p1_pdb, p2_pdb) pair; return (TMscore_P1, TMscore_P2)."""
    key = (os.path.abspath(p1_pdb), os.path.abspath(p2_pdb))
    if key not in cache:
        process = Popen(
            ["/install/USalign/USalign", p1_pdb, p2_pdb, '-mm', '0', '-ter', '2'],
            stdout=PIPE, stderr=PIPE,
        )
        stdout, _ = process.communicate()
        tmscore_list = [
            float(x.split(' ')[1])
            for x in stdout.decode().splitlines()
            if 'TM-score=' in x
        ]
        cache[key] = (
            tmscore_list[0] if tmscore_list else 0.0,
            tmscore_list[1] if len(tmscore_list) > 1 else 0.0,
        )
    return cache[key]


def main(args):
    P2 = args.target_pdb
    local_tmp_dir = os.getenv('TMPDIR')
    if not local_tmp_dir or not os.path.isdir(local_tmp_dir):
        local_tmp_dir = tempfile.mkdtemp(prefix='masif_mimicry_')
    else:
        os.makedirs(local_tmp_dir, exist_ok=True)
    
    if args.target_ppi_id not in ['p1', 'p2']:
        print('Please specify target_ppi_id as either p1 or p2')
        sys.exit(1)
    if args.target_ppi_id == 'p1':
        assert args.target_chain == P2.split('_')[1], f'Target chain {args.target_chain} does not match PDB {P2}'
    elif args.target_ppi_id == 'p2':
        assert args.target_chain == P2.split('_')[2], f'Target chain {args.target_chain} does not match PDB {P2}'

    if args.split_seed_list:
        lines = open(args.split_seed_list, 'r').readlines()
        lines = [line.strip() for line in lines]
        print(f'Number of proteins to search from: {len(lines)}')
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
    P2_all_feats = get_features(params, P2, args.target_ppi_id, source=False, flip_desc=False)

    if args.target_residue and args.target_chain and args.target_atom:
        try:
            P2_selected_points_idx = select_target_sites_by_radius_fps(
                np.array(P2_all_feats['mesh'].vertices),
                P2_all_feats['pdb'],
                chain=args.target_chain,
                residue=args.target_residue,
                atom_name=args.target_atom,
                radius=args.target_sampling_radius,
                num_points=args.num_points,
            )
        except ValueError as e:
            print(f'Error: {e}', flush=True)
            sys.exit(1)
        P2_patch_descs = P2_all_feats['desc'][P2_selected_points_idx]
        print(
            f'Searching sites similar to {P2} from {",".join(set(lines))}. \n'
            f'This will go through {len(P2_selected_points_idx)} points within '
            f'{args.target_sampling_radius} A of {args.target_atom} on residue '
            f'{args.target_residue} in {P2} chain {args.target_chain}...',
            flush=True,
        )
    else:
        P2_selected_points_idx, P2_patch_descs, P2_patch_iface = select_patches(
            P2_all_feats,
            downsample_rate=args.downsample,
            iface_cutoff=args.iface_cutoff,
            top_iface_percent=args.top_iface_percent,
            interface_only=args.interface_only,
            verbose=False
        )
        if len(P2_selected_points_idx) == 0:
            print(f'No points selected for {P2}. Please check the parameters. Exiting...')
            sys.exit(1)
        else:
            print(f'Searching sites similar to {P2} from {",".join(set(lines))}. \nWARNING: Exhausitive alignment mode. This will go through {len(P2_selected_points_idx)} points and will take some time...')

    scores = {}
    p2_output_root = os.path.join(args.output_dir, f'{P2}_{args.output_postfix}')
    os.makedirs(p2_output_root, exist_ok=True)
    write_target_vert_files(P2_all_feats, P2_selected_points_idx, args.target_ppi_id, p2_output_root)
    shutil.copy(P2_all_feats['pdb'], os.path.join(p2_output_root, f'{P2.split("_")[0]}_{args.target_chain}.pdb'))

    for P1 in lines:
        scores[(P1, P2)] = {
            'P1_id': [],
            'P2_id': [],
            'P1_source_ppi_id': [],
            'P2_source_ppi_id': [],
            'P1_source_site': [],
            'P2_source_site': [],
            'P1_source_TMscore': [],
            'P2_source_TMscore': [],
            'P1_source_iface': [],
            'P2_source_iface': [],
            'MaSIF-score': [],
            'flattened_transform': [],
        }
        if args.compute_source_residues:
            scores[(P1, P2)]['P1_source_residue'] = []
            scores[(P1, P2)]['P2_source_residue'] = []
        if args.count_clashes:
            scores[(P1, P2)]['ca_clash'] = []
            scores[(P1, P2)]['heavy_atom_clash'] = []

    if not args.interface_only:
        print(f'Aligning to all points with iface > {args.iface_cutoff} and desc dist < {args.desc_dist_cutoff}.')
    else:
        print(f'Aligning to interface points only with desc dist < {args.desc_dist_cutoff}.')

    tmscore_cache = {}
    target_structure_cache = {}

    for P1 in lines:

        if len(P1.split('_')) == 2:
            ppi_id_list = ['p1']
        elif len(P1.split('_')) == 3:
            ppi_id_list = ['p1', 'p2']
        else:
            raise ValueError('Invalid format. Please use PDB_X or PDB_X_X')

        for ppi_id in ppi_id_list:
            # try:
                P1_all_feats = get_features(params, P1, ppi_id, source=True, flip_desc=False)
                P1_selected_points_idx, P1_patch_descs, P1_patch_iface = select_patches(
                    P1_all_feats,
                    downsample_rate=args.downsample,
                    iface_cutoff=args.iface_cutoff,
                    top_iface_percent=args.top_iface_percent,
                    interface_only=args.interface_only,
                    verbose=False
                )

                # NOTE: compute a distance matrix between the selected points of P1 and P2 and filter them based on desc_dist_cutoff
                desc_dist = np.linalg.norm(P1_patch_descs[:, None, :] - P2_patch_descs[None, :, :], axis=-1)

                if np.min(desc_dist) > args.desc_dist_cutoff:
                    log(f'Skipping {P1} ({ppi_id}): no descriptors within {args.desc_dist_cutoff}')
                    continue

                seed_output_dir = os.path.join(p2_output_root, P1)
                log(f'Aligning {P1} ({ppi_id}) to {P2}...')
                total_matches = 0
                sites_aligned = 0

                if args.count_clashes:
                    P2_raw_pdb = os.path.join(
                        params['masif_target_root'],
                        'data_preparation', '00-raw_pdbs', f'{P2.split("_")[0]}.pdb',
                    )
                    P2_partner_chain_ids = parse_partner_chain_ids(P2, args.target_ppi_id)
                    P2_partner_chain_suffix = ''.join(P2_partner_chain_ids)
                    filtered_target_structure = get_filtered_target_structure(
                        P2_raw_pdb, P2_partner_chain_ids, target_structure_cache,
                    )
                else:
                    P2_raw_pdb = None
                    P2_partner_chain_suffix = None
                    filtered_target_structure = None

                for i, P2_center in enumerate(P2_selected_points_idx):
                    P1_selected_points_idx_final = P1_selected_points_idx[np.where(desc_dist[:, i] < args.desc_dist_cutoff)]

                    if len(P1_selected_points_idx_final) == 0:
                        log(f'Skipped site {P2_center} (no descriptors within {args.desc_dist_cutoff})')
                        continue

                    sites_aligned += 1
                    n_hits = 0

                    all_results, _, _, _ = multidock(
                        source_pt=P1_selected_points_idx_final,
                        source_pcd=P1_all_feats['pcd'], source_patch_idxs=P1_all_feats['indices'], source_descs=P1_all_feats['desc'], 
                        target_pt=P2_center,
                        target_pcd=P2_all_feats['pcd'], target_patch_idxs=P2_all_feats['indices'], target_descs=P2_all_feats['desc'], binder_align=False
                    )

                    output_root = os.path.join(args.output_dir, f'{P2}_{args.output_postfix}/{P1}')
                    for j, (result, P1_center) in enumerate(zip(all_results, P1_selected_points_idx_final)):
                        out_filename_base = f'{P1}_{ppi_id}_{P1_center}_to_{P2}_{args.target_ppi_id}_{P2_center}'

                        P1_patch_coords = transform_patch_coords(
                            P1_all_feats['pcd'],
                            P1_all_feats['indices'],
                            P1_center,
                            result.transformation,
                        )
                        descriptor_score = compute_descriptor_score(
                            P1_patch_coords,
                            P1_all_feats['desc'],
                            P1_all_feats['indices'],
                            P1_center,
                            P2_all_feats['pcd'],
                            P2_all_feats['desc'],
                            P2_all_feats['indices'],
                            P2_center,
                        )

                        if descriptor_score < args.desc_dist_score_cutoff:
                            continue

                        tmp_pdb_path = os.path.join(local_tmp_dir, f'{out_filename_base}.pdb')
                        _ = transform_structure(
                            P1_all_feats['pdb'], result.transformation, tmp_pdb_path,
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
                                if os.path.exists(tmp_pdb_path):
                                    os.remove(tmp_pdb_path)
                                continue
                        else:
                            ca_clashes, heavy_clashes = 0, 0
                            masif_score = descriptor_score

                        n_hits += 1
                        os.makedirs(output_root, exist_ok=True)
                        out_pdb_path = os.path.join(output_root, f'{out_filename_base}.pdb')
                        shutil.move(tmp_pdb_path, out_pdb_path)

                        if args.count_clashes:
                            io = PDBIO()
                            io.set_structure(filtered_target_structure)
                            io.save(os.path.join(
                                args.output_dir, f'{P2}_{args.output_postfix}',
                                f'{P2.split("_")[0]}_{P2_partner_chain_suffix}.pdb',
                            ))

                        if args.compute_source_residues:
                            target_atom, _ = surf2atom(
                                point_coords=np.array(P2_all_feats['pcd'].points)[P2_center].reshape(1, -1),
                                pdb_path=P2_all_feats['pdb'],
                            )
                            target_residue = target_atom[0].get_parent().get_id()[1]
                            P1_center_coord = apply_transform(
                                np.asarray(P1_all_feats['pcd'].points)[P1_center], result.transformation,
                            )
                            P1_nearest_atom, _ = surf2atom(
                                point_coords=P1_center_coord.reshape(1, -1),
                                pdb_path=out_pdb_path,
                            )
                            P1_nearest_res = P1_nearest_atom[0].get_parent().get_id()[1]

                        TMscore_P1, TMscore_P2 = get_usalign_tmscores(
                            P1_all_feats['pdb'], P2_all_feats['pdb'], tmscore_cache,
                        )

                        scores[(P1, P2)]['P1_id'].append(P1)
                        scores[(P1, P2)]['P2_id'].append(P2)
                        scores[(P1, P2)]['P1_source_ppi_id'].append(ppi_id)
                        scores[(P1, P2)]['P2_source_ppi_id'].append(args.target_ppi_id)
                        if args.compute_source_residues:
                            scores[(P1, P2)]['P1_source_residue'].append(P1_nearest_res)
                            scores[(P1, P2)]['P2_source_residue'].append(target_residue)
                        scores[(P1, P2)]['P1_source_site'].append(P1_center)
                        scores[(P1, P2)]['P2_source_site'].append(P2_center)
                        scores[(P1, P2)]['P1_source_TMscore'].append(TMscore_P1)
                        scores[(P1, P2)]['P2_source_TMscore'].append(TMscore_P2)
                        scores[(P1, P2)]['P1_source_iface'].append(P1_all_feats['iface'][0][P1_center])
                        scores[(P1, P2)]['P2_source_iface'].append(P2_all_feats['iface'][0][P2_center])
                        scores[(P1, P2)]['MaSIF-score'].append(masif_score)
                        scores[(P1, P2)]['flattened_transform'].append(flatten_transform(result.transformation))
                        if args.count_clashes:
                            scores[(P1, P2)]['ca_clash'].append(ca_clashes)
                            scores[(P1, P2)]['heavy_atom_clash'].append(heavy_clashes)

                    total_matches += n_hits
                    log(f'Found {n_hits} matches to site {P2_center}')

                log(f'Done {P1} ({ppi_id}): {total_matches} matches across {sites_aligned} sites')

                if len(scores[(P1, P2)]['P1_id']) > 0:
                    pd.DataFrame(scores[(P1, P2)]).to_csv(
                        f'{seed_output_dir}/{P1}_{ppi_id}_to_{P2}_{args.target_ppi_id}.csv',
                        index=False,
                    )
            # except Exception as e:
            #     print(f'Error: {e}')

if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    main(args)
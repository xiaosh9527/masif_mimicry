import os, sys, glob, shutil, pickle, csv, argparse, copy
import numpy as np
import pandas as pd

from utils import *
from subprocess import Popen, PIPE

def create_parser():
    parser = argparse.ArgumentParser(description='MaSIF Pipeline for mimicry search.')
    parser.add_argument('--split_seed_list', type=str, help='List that specify uniprot ID of the proteins to search from.')
    parser.add_argument('--seed_db_root', type=str, default='/work/lpdi/users/shxiao/masif_seed/masif/', help='Root directory where seed PDB is processed')
    parser.add_argument('--seed_db', type=str, help='database name where seed PDB is processed')
    parser.add_argument('--seed_pdb', type=str, help='specify uniprot ID of the proteins to search from.')
    parser.add_argument('--target_db_root', type=str, default='/work/lpdi/users/shxiao/masif_seed/masif_seed_search', help='Root directory where target PDB is processed')
    parser.add_argument('--target_db', type=str, default='masif_degron', help='database name where target PDB is processed')
    parser.add_argument('--target_pdb', type=str, help='Target PDB ID')
    parser.add_argument('--target_ppi_id', type=str, default='p1', help='Target PPI ID')
    parser.add_argument('--target_chain', type=str, help='Target PDB Chain ID')
    parser.add_argument('--target_residue', type=int, help='Target residue ID')
    parser.add_argument('--target_atom', type=str, help='Target atom name')
    parser.add_argument('--top_iface_percent', type=float, default=0.02, help='Top interface percent to select points from')
    parser.add_argument('--downsample', type=int, default=1, help='Downsample rate for the sites on proteome proteins')
    parser.add_argument('--output_dir', type=str, default='.', help='Output directory')
    parser.add_argument('--output_postfix', type=str, default='outputs', help='Output postfix')
    parser.add_argument('--num_points', type=int, default=5, help='Number of nearest points to search')
    parser.add_argument('--masif_app', type=str, default='ppi_search', help='MaSIF applications: can do ppi_search or site')
    parser.add_argument('--iface_cutoff', type=float, default=0.0, help='Interface cutoff')
    parser.add_argument('--interface_only', action='store_true', help='Only align to interfaces in the seed database; this requires the seed database to be protein complexes.')
    parser.add_argument('--desc_dist_cutoff', type=float, default=1.5, help='Desc dist cutoff')
    parser.add_argument('--desc_dist_score_cutoff', type=float, default=0.25, help='Desc dist score cutoff. Recommend 30 for masif_app ppi_search and 15 for site.')
    parser.add_argument('--count_clashes', action='store_true', help='Count clashes')
    parser.add_argument('--ca_clash_threshold', type=int, default=1, help='CA clash threshold')
    parser.add_argument('--heavy_atom_clash_threshold', type=int, default=5, help='Clash threshold')
    return parser

def main(args):
    P2 = args.target_pdb
    local_tmp_dir = os.getenv('TMPDIR')    
    
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

    print(f'Searching sites similar to {P2} from {",".join(set(lines))}...')
        
    params = set_params(target_root = f'{args.target_db_root}/data/{args.target_db}', masif_db_root = args.seed_db_root, db_name = args.seed_db, masif_app = args.masif_app)
    P2_all_feats = get_features(params, P2, args.target_ppi_id, source=False, flip_desc=False)

    if args.target_residue and args.target_chain and args.target_atom:
        P2_patch_iface = P2_all_feats['iface'][0]
        P2_nearest_points_idx = res2surf(np.array(P2_all_feats['mesh'].vertices), P2_all_feats['pdb'], chain = args.target_chain, residue = args.target_residue, atom_name = args.target_atom, k = 10)
        P2_nearest_points_iface = P2_patch_iface[P2_nearest_points_idx[0]]
        P2_selected_points_idx = P2_nearest_points_idx[0][P2_nearest_points_iface.argsort()[::-1][:args.num_points]]
        P2_patch_descs = []
        for idx in P2_selected_points_idx:
            P2_patch_descs.append(P2_all_feats['desc'][idx])
        P2_patch_descs = np.array(P2_patch_descs)
        print(f'This will go through {len(P2_selected_points_idx)} points near {args.target_residue} in {P2} chain {args.target_chain}...')
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
            print(f'WARNING: Exhausitive alignment mode. This will go through {len(P2_selected_points_idx)} points and will take some time...')

    scores = {}
    os.makedirs(os.path.join(args.output_dir, f'{P2}_{args.output_postfix}'), exist_ok=True)
    shutil.copy(P2_all_feats['pdb'], os.path.join(args.output_dir, f'{P2}_{args.output_postfix}', f'{P2}.pdb'))

    for P1 in lines:
        scores[(P1, P2)] = {
            'P1_id': [],
            'P2_id': [],
            'P1_source_ppi_id': [],
            'P2_source_ppi_id': [],
            'P1_source_residue': [],
            'P2_source_residue': [],
            'P1_source_site': [],
            'P2_source_site': [],
            'P1_source_TMscore': [],
            'P2_source_TMscore': [],
            'P1_source_iface': [],
            'P2_source_iface': [],
            'MaSIF-score': [],
        }
        if args.count_clashes:
            scores[(P1, P2)]['ca_clash'] = []
            scores[(P1, P2)]['heavy_atom_clash'] = []

    if not args.interface_only:
        print(f'Aligning to all points with iface > {args.iface_cutoff} and desc dist < {args.desc_dist_cutoff}...')
    else:
        print(f'Aligning to interface points only with desc dist < {args.desc_dist_cutoff}...')

    for P1 in lines:
        if len(P1.split('_')) == 2:
            ppi_id_list = ['p1']
        elif len(P1.split('_')) == 3:
            ppi_id_list = ['p1', 'p2']
        else:
            raise ValueError('Invalid format. Please use PDB_X or PDB_X_X')

        for ppi_id in ppi_id_list:
            try:
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

                if np.min(desc_dist) > args.desc_dist_cutoff: continue

                for i, P2_center in enumerate(P2_selected_points_idx):  
                    P1_selected_points_idx_final = P1_selected_points_idx[np.where(desc_dist[:, i]<args.desc_dist_cutoff)]

                    if len(P1_selected_points_idx_final) == 0: continue
                    
                    target_atom, _ = surf2atom(
                        point_coords = np.array(P2_all_feats['pcd'].points)[P2_center].reshape(1,-1),
                        pdb_path = P2_all_feats['pdb'],
                        k=1,
                        distance_cutoff=5.0, 
                        rsa_cutoff = None, 
                        exclude_backbone = True
                    )
                    
                    target_residue = target_atom[0].get_parent().get_id()[1]

                    all_results, _, _, _ = multidock(
                        source_pt=P1_selected_points_idx_final,
                        source_pcd=P1_all_feats['pcd'], source_patch_idxs=P1_all_feats['indices'], source_descs=P1_all_feats['desc'], 
                        target_pt=P2_center,
                        target_pcd=P2_all_feats['pcd'], target_patch_idxs=P2_all_feats['indices'], target_descs=P2_all_feats['desc'], binder_align=False
                    )

                    for j, (result, P1_center) in enumerate(zip(all_results, P1_selected_points_idx_final)):
                        output_root = os.path.join(args.output_dir, f'{P2}_{args.output_postfix}/{P1}')
                        out_filename_base = f'{P1}_{ppi_id}_{P1_center}_to_{P2}_{args.target_ppi_id}_{P2_center}'

                        P1_pcd = copy.deepcopy(P1_all_feats['pcd'])
                        P1_pcd.transform(result.transformation)
                        _ = transform_structure(P1_all_feats['pdb'], result.transformation, os.path.join(local_tmp_dir, f'{out_filename_base}.pdb'))

                        if args.count_clashes:
                            P2_raw_pdb = os.path.join(params['masif_target_root'], 'data_preparation', '00-raw_pdbs', f'{P2.split("_")[0]}.pdb')
                            P2_p2_chain = P2.split('_')[-1]
                            output = compute_score_and_clashes(
                                P1_pdb = os.path.join(local_tmp_dir, f'{out_filename_base}.pdb'), P1_pcd = P1_pcd, P1_descs = P1_all_feats['desc'], P1_site = P1_center, P1_indices = P1_all_feats['indices'],
                                P2_pdb = P2_all_feats['pdb'], P2_pcd = P2_all_feats['pcd'], P2_descs = P2_all_feats['desc'], P2_site = P2_center, P2_indices = P2_all_feats['indices'],
                                compute_clashes = True, target_structure = P2_raw_pdb, target_chain = P2_p2_chain, ca_clash_threshold = args.ca_clash_threshold, heavy_atom_clash_threshold = args.heavy_atom_clash_threshold
                            )
                        else:
                            output = compute_score_and_clashes(
                                P1_pdb = os.path.join(local_tmp_dir, f'{out_filename_base}.pdb'), P1_pcd = P1_pcd, P1_descs = P1_all_feats['desc'], P1_site = P1_center, P1_indices = P1_all_feats['indices'],
                                P2_pdb = P2_all_feats['pdb'], P2_pcd = P2_all_feats['pcd'], P2_descs = P2_all_feats['desc'], P2_site = P2_center, P2_indices = P2_all_feats['indices'],
                                compute_clashes = False
                            )

                        if output[1] >= args.desc_dist_score_cutoff and output[0][0] <= args.ca_clash_threshold and output[0][1] <= args.heavy_atom_clash_threshold:                    
                            os.makedirs(output_root, exist_ok=True)
                            shutil.move(os.path.join(local_tmp_dir, f'{out_filename_base}.pdb'), os.path.join(output_root, f'{out_filename_base}.pdb'))
                            P1_raw_pdb = os.path.join(params['top_seed_dir'], 'data_preparation', '00-raw_pdbs', f'{P1.split("_")[0]}.pdb')
                            _ = transform_structure(P1_raw_pdb, result.transformation, os.path.join(output_root, f'{P1.split("_")[0]}.pdb'))
                            
                            # NOTE: Find the nearest residue on P1 to the center point
                            P1_nearest_atom, _ = surf2atom(
                                point_coords = np.array(P1_pcd.points)[P1_center].reshape(1,-1), 
                                pdb_path = os.path.join(output_root, f'{out_filename_base}.pdb'), 
                                k = 1, 
                                distance_cutoff = 3.0, 
                                rsa_cutoff = 0.2, 
                                exclude_backbone = True
                            )

                            P1_nearest_res = P1_nearest_atom[0].get_parent().get_id()[1]

                            # Compute TMscore using USalign
                            process = Popen(["/work/lpdi/users/shxiao/USalign", P1_all_feats['pdb'], P2_all_feats['pdb'], '-mm', '0', '-ter', '2'], stdout=PIPE, stderr=PIPE)
                            stdout, _ = process.communicate()
                            TMscore_P1 = float(stdout.decode().splitlines()[14].split(' ')[1])
                            TMscore_P2 = float(stdout.decode().splitlines()[15].split(' ')[1])

                            scores[(P1, P2)]['P1_id'].append(P1)
                            scores[(P1, P2)]['P2_id'].append(P2)
                            scores[(P1, P2)]['P1_source_ppi_id'].append(ppi_id)
                            scores[(P1, P2)]['P2_source_ppi_id'].append(args.target_ppi_id)
                            scores[(P1, P2)]['P1_source_residue'].append(P1_nearest_res)
                            scores[(P1, P2)]['P2_source_residue'].append(target_residue)
                            scores[(P1, P2)]['P1_source_site'].append(P1_center)
                            scores[(P1, P2)]['P2_source_site'].append(P2_center)
                            scores[(P1, P2)]['P1_source_TMscore'].append(TMscore_P1)
                            scores[(P1, P2)]['P2_source_TMscore'].append(TMscore_P2)
                            scores[(P1, P2)]['P1_source_iface'].append(P1_all_feats['iface'][0][P1_center])
                            scores[(P1, P2)]['P2_source_iface'].append(P2_all_feats['iface'][0][P2_center])
                            scores[(P1, P2)]['MaSIF-score'].append(output[1])
                            if args.count_clashes:
                                scores[(P1, P2)]['ca_clash'].append(output[0][0])
                                scores[(P1, P2)]['heavy_atom_clash'].append(output[0][1])

                    if len(scores[(P1, P2)]['P1_id'])>0:
                        pd.DataFrame(scores[(P1, P2)]).to_csv(f'{output_root}/{P1}_{ppi_id}_to_{P2}_{args.target_ppi_id}.csv')
                    else:
                        pass

            except Exception as e:
                print(f'Error: {e}')

if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    main(args)
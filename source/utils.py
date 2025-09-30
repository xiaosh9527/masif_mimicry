import os
import numpy as np
import open3d as o3d
from scipy.spatial import cKDTree
from Bio.PDB import PDBParser, PDBIO, Selection, DSSP

from geometry.open3d_import import *

from default_config.masif_opts import masif_opts

def set_params(target_root, masif_db_root: str = "/work/lpdi/users/shxiao/scratch/masif_seed/masif", db_name: str = "masif_oas", masif_app: str = "ppi_search"):
    params = {}
    # Directory where the database is located.
    params["masif_db_root"] = masif_db_root
    # Seeds (i.e., the fragments) that will be used for this search.
    params["top_seed_dir"] = os.path.join(params["masif_db_root"], f"data/{db_name}/")
    # Root of the targets directory (where to find the sources.)
    params["masif_target_root"] = target_root
    # Output directory (target_name, target_site, matched_seed)
    params["out_dir_template"] = "tmp/{}/"

    # Seed locations
    params["seed_surf_dir"] = os.path.join(params["top_seed_dir"], masif_opts["ply_chain_dir"])
    params["seed_iface_dir"] = os.path.join(params["top_seed_dir"], masif_opts["site"]["out_pred_dir"])
    params["seed_ply_iface_dir"] = os.path.join(params["top_seed_dir"], masif_opts["site"]["out_surf_dir"])
    params["seed_pdb_dir"] = os.path.join(params["top_seed_dir"], masif_opts["pdb_chain_dir"])
    params["seed_desc_dir"] = os.path.join(params["top_seed_dir"], masif_opts["ppi_search"]["desc_dir"])
    params["seed_precomp_dir"] = os.path.join(params["top_seed_dir"], masif_opts[masif_app]["masif_precomputation_dir"])

    # Target locations
    params["target_surf_dir"] = os.path.join(params["masif_target_root"], masif_opts["ply_chain_dir"])
    params["target_iface_dir"] = os.path.join(params["masif_target_root"], masif_opts["site"]["out_pred_dir"])
    params["target_ply_iface_dir"] = os.path.join(params["masif_target_root"], masif_opts["site"]["out_surf_dir"])
    params["target_pdb_dir"] = os.path.join(params["masif_target_root"], masif_opts["pdb_chain_dir"])
    params["target_desc_dir"] = os.path.join(params["masif_target_root"], masif_opts["ppi_search"]["desc_dir"])
    params["target_precomp_dir"] = os.path.join(params["masif_target_root"], masif_opts[masif_app]["masif_precomputation_dir"])

    return params

def get_features(params, pdb: str, pid: str, source: bool = True, flip_desc: bool = False):
    '''
    Load all the features for a given pdb and chain id.
    pdb: PDB code with chain info (e.g., 1a2k_A or 1a2k_A_B)
    pid: "p1" or "p2"
    source: Whether this is the source (seed) or target (receptor)
    flip_desc: Whether to use the descriptors with flipped input features. 
    Returns a dictionary with all the features.
    '''
    
    if pid == "p1":
        P = pdb.split("_")[0] + "_" + pdb.split("_")[1]
    else:
        P = pdb.split("_")[0] + "_" + pdb.split("_")[2]

    if source:
        mask = "seed"
    else:
        mask = "target"

    pdb_fn = os.path.join(params[f"{mask}_pdb_dir"], f"{P}.pdb")
    ply = os.path.join(params[f"{mask}_surf_dir"], f"{P}.ply")
    input_feat_fn = os.path.join(params[f"{mask}_precomp_dir"], pdb, f"{pid}_input_feat.npy")
    indices_fn = os.path.join(params[f"{mask}_precomp_dir"], pdb, f"{pid}_list_indices.npy")
    rho_fn = os.path.join(params[f"{mask}_precomp_dir"], pdb, f"{pid}_rho_wrt_center.npy")
    theta_fn = os.path.join(params[f"{mask}_precomp_dir"], pdb, f"{pid}_theta_wrt_center.npy")
    ilabel_fn = os.path.join(params[f"{mask}_precomp_dir"], pdb, f"{pid}_iface_labels.npy")
    iface_fn = os.path.join(params[f"{mask}_iface_dir"], "pred_" + P + ".npy")
    if flip_desc:
        desc_fn = os.path.join(params[f"{mask}_desc_dir"], pdb, f"{pid}_desc_flipped.npy")
    else:
        desc_fn = os.path.join(params[f"{mask}_desc_dir"], pdb, f"{pid}_desc_straight.npy")
    
    all_feats = dict(
        pdb = pdb_fn,
        mesh = o3d.io.read_triangle_mesh(ply),
        pcd = o3d.io.read_point_cloud(ply),
        rho = np.load(rho_fn),
        theta = np.load(theta_fn),
        desc = np.load(desc_fn),
        input_feat = np.load(input_feat_fn),
        indices = np.load(indices_fn, allow_pickle=True),
        iface = np.load(iface_fn),
        ilabel = np.load(ilabel_fn),
    )

    return all_feats

def surf2atom(
    point_coords,
    pdb_path: str,
    k: int = 1,
    distance_cutoff: float = 5.0, 
    rsa_cutoff: float = 0.1, 
    exclude_backbone: bool = True,
    dssp_bin_path: str = "/work/lpdi/bin/sequence/dssp"
    ):
    '''
    This function takes the pdb file and the surface point coordinates and returns the closest surface-exposed atoms.
    point_coords: Nx3 numpy array of surface point coordinates
    pdb_path: Path to the PDB file
    k: Number of closest atoms to return
    distance_cutoff: Maximum distance to consider for closest atoms
    rsa_cutoff: Minimum relative solvent accessibility to consider an atom as surface-exposed
    exclude_backbone: Whether to exclude backbone atoms (N, CA, C, O)
    dssp_bin_path: Path to the DSSP binary

    Returns:
        nearest_atom: List of closest atoms to each point (list of Bio.PDB.Atom objects)
        nearest_atom_ss: List of secondary structure labels of the closest atoms (H, E, C)
    '''
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('', pdb_path)

    # NOTE: Defining surface atoms as RSA higher than 0.1
    dssp = DSSP(struct[0], pdb_path, dssp=dssp_bin_path)
    ss_atom_list = np.concatenate([[x[2]]*len(list(y.get_atoms())) for x,y in zip(dssp, struct.get_residues())])
    rsa_atom_list = np.concatenate([[x[3]]*len(list(y.get_atoms())) for x,y in zip(dssp, struct.get_residues())])
    rsa_atom_list = np.array([float(x) if x != 'NA' else 1.0 for x in rsa_atom_list])

    if rsa_cutoff is not None:
        if exclude_backbone:
            atoms = np.array([a for a,b in zip(struct.get_atoms(),rsa_atom_list) if b >= rsa_cutoff and a.get_id() not in ['N', 'CA', 'C', 'O']])
            ss_label = np.array([c for a,b,c in zip(struct.get_atoms(),rsa_atom_list,ss_atom_list) if b >= rsa_cutoff and a.get_id() not in ['N', 'CA', 'C', 'O']])
        else:
            atoms = np.array([a for a,b in zip(struct.get_atoms(),rsa_atom_list) if b >= rsa_cutoff])
            ss_label = np.array([a for a,b in zip(ss_atom_list,rsa_atom_list) if b >= rsa_cutoff])
    else:
        if exclude_backbone:
            atoms = np.array([a for a in struct.get_atoms() if a.get_id() not in ['N', 'CA', 'C', 'O']])
            ss_label = np.array([a for a,b in zip(ss_atom_list, struct.get_atoms()) if b.get_id() not in ['N', 'CA', 'C', 'O']])
        else:
            atoms = np.array(list(struct.get_atoms()))
            ss_label = ss_atom_list

    atom_coords = np.array([atom.get_coord() for atom in atoms])

    # NOTE: this will register the points to the closest atoms/residues
    atom_cdktree = cKDTree(atom_coords)
    d, nearest_atoms_idx = atom_cdktree.query(point_coords, k=k)
        
    if k == 1:
        nearest_atom = [a for a in atoms[nearest_atoms_idx]]
        nearest_atom_ss = [a for a in ss_label[nearest_atoms_idx]]
    else:
        nearest_atom = []
        nearest_atom_ss = []
        for i in range(d.shape[0]):
            nearest_atom.append([atoms[nearest_atoms_idx[i, j]] for j in range(d.shape[1]) if d[i, j] <= distance_cutoff])
            nearest_atom_ss.append([ss_label[nearest_atoms_idx[i, j]] for j in range(d.shape[1]) if d[i, j] <= distance_cutoff])
            
    return nearest_atom, nearest_atom_ss
    
def res2surf(point_coords, pdb_path: str, chain: str, residue: int, atom_name: str, k: int = 1):
    '''
    This function takes the pdb file and the surface point coordinates and returns the closest surface points to a given residue/atom
    
    point_coords: Nx3 numpy array of surface point coordinates
    pdb_path: Path to the PDB file
    chain: Chain identifier of the residue
    residue: Residue number
    atom_name: Atom name (e.g., 'CA', 'CB', 'N',
                'C', 'O', 'CG', etc.)   
    k: Number of closest points to return

    Returns:
        nearest_points_idx: Indices of the closest surface points to the specified atom
    '''
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('target', pdb_path)
    atoms = [atom for atom in struct.get_atoms() if atom.get_parent().get_id()[1]==residue and atom.get_parent().get_parent().get_id()==chain and atom.get_id()==atom_name]
    atom_coords = np.array([atom.get_coord() for atom in atoms])
    
    # NOTE: this will register the atom/residue to the closest points
    point_cdktree = cKDTree(point_coords)
    d, nearest_points_idx = point_cdktree.query(atom_coords, k=k)
    
    return nearest_points_idx

def select_patches(
    P_all_feats, 
    downsample_rate: int = 5, 
    iface_cutoff: float = 0.0, 
    min_patch_num: int = 50, 
    top_iface_percent: float = None, 
    interface_only: bool = False,
    verbose: bool = True
):  
    """
    Select patches from the point cloud based on indices and downsample rate.
    P_all_feats: Dictionary containing all features of the point cloud.
    downsample_rate: Downsample rate for selecting points from each patch.
    iface_cutoff: Minimum interface score to consider a patch.
    min_patch_num: Minimum number of patches to select.
    top_iface_percent: If specified, select only the top x percentage of patches based on interface score.
    interface_only: If True, select only patches that are labeled as interface.
    verbose: Whether to print information about the selection process.
    
    Returns:
        selected_points_idx: Indices of the selected patches.
        patch_descs: Descriptors of the selected patches.
        patch_iface: Interface scores of the selected patches.
    """
    subpcd_coverage = set()
    selected_points_idx = []
    patch_descs = []
    patch_iface = P_all_feats['iface'][0]
    P_indices = P_all_feats['indices']
    P_interface_points = np.where(P_all_feats['ilabel'] == 1)[0]
    for ii in range(len(P_all_feats['desc'])):
        desc = P_all_feats['desc'][ii]
        if interface_only:
            if ii not in subpcd_coverage and ii in P_interface_points:
                subpcd_coverage.update(P_indices[ii][:downsample_rate])
                selected_points_idx.append(ii)
                patch_descs.append(desc)
        else:
            if ii not in subpcd_coverage and patch_iface[ii] >= iface_cutoff:
                subpcd_coverage.update(P_indices[ii][:downsample_rate])
                selected_points_idx.append(ii)
                patch_descs.append(desc)
    
    patch_descs = np.array(patch_descs)
    selected_points_idx = np.array(selected_points_idx)
    
    if len(selected_points_idx) == 0:
        print("No points selected. Please check the parameters.")
        return [], [], []
    
    if top_iface_percent is not None:
        top_iface_num = round(top_iface_percent*len(P_all_feats['desc']))
        top_iface_num = max(top_iface_num, min_patch_num)
        top_iface_idx = np.argsort(patch_iface[selected_points_idx])[::-1][:top_iface_num]
    else:
        top_iface_idx = np.arange(len(selected_points_idx))
        
    selected_points_idx = selected_points_idx[top_iface_idx]
    patch_descs = patch_descs[top_iface_idx]

    if verbose:
        if interface_only:
            print(f'WARNING: Exhausitive alignment mode. This will go through {len(selected_points_idx)} interface points and will take some time...')
        else:
            print(f'Selected {len(selected_points_idx)} points with interface score >= {iface_cutoff} and downsample rate {downsample_rate}.')
        
    return selected_points_idx, patch_descs, patch_iface[selected_points_idx]

def get_patch_geo(
        pcd,
        patch_coords,
        center,
        descriptors,
        outward_shift=0.25,
        flip_normals=False):
    """
    Returns a patch from a point cloud pcd with center point center (int),
    based on geodesic distances from patch coords and corresponding Feature descriptors.
    """
    patch_idxs = patch_coords[center]
    patch_pts = np.asarray(pcd.points)[patch_idxs, :]
    patch_nrmls = np.asarray(pcd.normals)[patch_idxs, :]
    patch_pts = patch_pts + outward_shift * patch_nrmls
    if flip_normals:
        patch_nrmls = -patch_nrmls

    patch = o3d.geometry.PointCloud()
    patch.points = o3d.utility.Vector3dVector(patch_pts)
    patch.normals = o3d.utility.Vector3dVector(patch_nrmls)
    patch_descs = [o3d.pipelines.registration.Feature(), o3d.pipelines.registration.Feature(), o3d.pipelines.registration.Feature()]
    patch_descs[0].data = descriptors[patch_idxs, :].T
    
    return patch, patch_descs, patch_idxs
    
def multidock(
    source_pt, source_pcd, source_patch_idxs, source_descs, 
    target_pt, target_pcd, target_patch_idxs, target_descs, 
    binder_align: bool = False, ransac_skip: bool = False
):
    ransac_radius=1.5
    ransac_iter=1000000000
    all_results = []
    all_source_patch = []
    all_source_scores = []
    all_source_desc = []
    all_source_idx = []
    
    target_patch, target_patch_descs, target_patch_idx = get_patch_geo(
            target_pcd, target_patch_idxs, target_pt, target_descs, flip_normals=binder_align, outward_shift=0.25)
    
    for pt in source_pt:
        source_patch, source_patch_descs, source_patch_idx = get_patch_geo(
            source_pcd, source_patch_idxs, pt, source_descs, outward_shift=0.25)
        
        if not ransac_skip: 
            result = o3d.pipelines.registration.registration_ransac_based_on_feature_matching(
                source=source_patch, target=target_patch, source_feature=source_patch_descs[0], target_feature=target_patch_descs[0], 
                mutual_filter=False, max_correspondence_distance=ransac_radius, 
                estimation_method=o3d.pipelines.registration.TransformationEstimationPointToPoint(False), 
                ransac_n = 3, checkers=[o3d.pipelines.registration.CorrespondenceCheckerBasedOnEdgeLength(0.9),
                o3d.pipelines.registration.CorrespondenceCheckerBasedOnDistance(1.0),
                o3d.pipelines.registration.CorrespondenceCheckerBasedOnNormal(np.pi/2)],
                criteria=o3d.pipelines.registration.RANSACConvergenceCriteria(max_iteration=ransac_iter, confidence=0.999999),
                seed=42
            )
            # result = registration_ransac_based_on_feature_matching(
            #     source=source_patch, target=target_patch, source_feature=source_patch_descs[0], target_feature=target_patch_descs[0], 
            #     max_correspondence_distance=ransac_radius,
            #     estimation_method=TransformationEstimationPointToPoint(False), ransac_n=3,
            #     checkers=[CorrespondenceCheckerBasedOnEdgeLength(0.9),
            #     CorrespondenceCheckerBasedOnDistance(1.0),
            #     CorrespondenceCheckerBasedOnNormal(np.pi/2)],
            #     criteria=RANSACConvergenceCriteria(ransac_iter, 500)
            # )
            init = result.transformation
        else:
            init = np.identity(4)

        result_icp = o3d.pipelines.registration.registration_icp(
            source=source_patch, target=target_patch,
            max_correspondence_distance=1.5, init=init,
            estimation_method=o3d.pipelines.registration.TransformationEstimationPointToPlane(),
            criteria=o3d.pipelines.registration.ICPConvergenceCriteria()
        )
        # result_icp = registration_icp(
        #     source=source_patch, target=target_patch,
        #     max_correspondence_distance=1.5, init=init, 
        #     estimation_method=TransformationEstimationPointToPlane(),
        #     criteria=ICPConvergenceCriteria()
        # )
                
        source_patch.transform(result_icp.transformation)
        all_results.append(result_icp)
        all_source_patch.append(source_patch)
        all_source_desc.append(source_patch_descs)
        all_source_idx.append(source_patch_idx)

    return all_results, all_source_patch, all_source_desc, all_source_idx 

def compute_score_and_clashes(
    P1_pdb, P1_pcd, P1_descs, P1_site, P1_indices, 
    P2_pdb, P2_pcd, P2_descs, P2_site, P2_indices, 
    compute_clashes: bool,
    **kwargs
):  

    '''
    Compute the descriptor distance score and the number of clashes between two structures.
    P1_pdb: Path to the PDB file of structure 1 (source).
    P1_pcd: Open3D point cloud of structure 1.
    P1_descs: NxD numpy array of descriptors for structure 1.
    P1_site: Index of the site on structure 1.
    P1_indices: List of indices mapping surface points to original point cloud for structure 1
    P2_pdb: Path to the PDB file of structure 2 (target).
    P2_pcd: Open3D point cloud of structure 2.
    P2_descs: NxD numpy array of descriptors for structure 2.
    P2_site: Index of the site on structure 2.
    P2_indices: List of indices mapping surface points to original point cloud for structure 2
    compute_clashes: Whether to compute clashes.
    kwargs: Additional arguments for clash computation (e.g., target_structure, target_chain).
    
    Returns:
        desc_dist_score: Descriptor distance score.
        clashing_ca: Number of clashing CA atoms (if compute_clashes is True).
        clashing: Number of clashing heavy atoms (if compute_clashes is True).
    '''
    
    P1_patch_coords = np.array(P1_pcd.points)[P1_indices[P1_site]]
    P2_patch_coords = np.array(P2_pcd.points)[P2_indices[P2_site]]

    P2_patch_ckdtree = cKDTree(P2_patch_coords)
    d_nn, r_nn = P2_patch_ckdtree.query(P1_patch_coords)
    neigh = np.where(d_nn <= 5.0)[0]

    if len(neigh) == 0:
        return [[np.inf, np.inf], 0.0]

    P1_patch_desc = P1_descs[P1_indices[P1_site]]
    P2_patch_desc = P2_descs[P2_indices[P2_site]]

    desc_dist_score = np.sum(np.square((1 / np.sqrt(np.sum(np.square(P1_patch_desc - P2_patch_desc[r_nn]), axis=1)))[neigh]))
    normalized_desc_dist_score = np.tanh(desc_dist_score/80)
    output = [[0, 0], normalized_desc_dist_score]
    
    if compute_clashes:
        pdb_parser = PDBParser(QUIET=True)
        assert 'target_structure' in kwargs, 'Please provide target_structure for clash counting.'
        assert 'target_chain' in kwargs, 'Please provide target_chain for clash counting.'
        source_structure = pdb_parser.get_structure('', P1_pdb) # NOTE: P1_pdb is already transformed
        target_structure = pdb_parser.get_structure('', kwargs['target_structure'])
        target_chain = list(kwargs['target_chain'])
        
        chains_to_remove = []
        for chain in target_structure.get_chains():
            if chain.id not in target_chain:
                chains_to_remove.append(chain)
        
        for chain in chains_to_remove:
            target_structure[0].detach_child(chain.id)

        output[0][0], output[0][1] = count_clashes(source_structure, target_structure, radius=2.0)
    return output

def count_clashes(
    source_structure,
    target_structure, 
    radius=2.0, 
    ):
    '''
    Count the number of clashing atoms between source and target structures.
    Clashing is defined as having any atom within a certain radius (default 2.0 A).
    Returns the number of clashing CA atoms and the number of clashing heavy atoms.

    source_structure: Bio.PDB.Structure object of the source (transformed) structure.
    target_structure: Bio.PDB.Structure object of the target structure.
    radius: distance threshold to consider a clash.

    Returns:
    clashing_ca: Number of clashing CA atoms.
    clashing: Number of clashing heavy atoms.
    '''

    source_ca_coords = np.array([atom.get_coord() for atom in source_structure.get_atoms() if atom.get_id() == 'CA'])

    target_ca_coords = np.array([atom.get_coord() for atom in target_structure.get_atoms() if atom.get_id() == 'CA'])
    target_pcd_tree = cKDTree(target_ca_coords)

    d_nn_ca, _ = target_pcd_tree.query(np.asarray(source_ca_coords), k=1, distance_upper_bound=radius)
    clashing_ca = np.sum(d_nn_ca<=radius)
    
    source_atoms = [atom for atom in source_structure.get_atoms() if not atom.get_name().startswith('H')]
    source_coords = np.array([atom.get_coord() for atom in source_atoms])

    d_nn, _ = target_pcd_tree.query(np.asarray(source_coords), k=1, distance_upper_bound=radius)
    clashing = np.sum(d_nn<=radius)
    
    return clashing_ca, clashing

def transform_structure(input_path, T, output_path=None):
    '''
    Transform the 3D structure by applying a transformation matrix T.
    
    input_path: Path to the input PDB file.
    T: 4x4 transformation matrix.
    output_path: Path to save the transformed PDB file. If None, the file is not saved.
    
    Returns the transformed structure as a Bio.PDB Structure object.
    '''
    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure('', input_path)
    P_atoms = [atom for atom in structure.get_atoms() if not atom.get_name().startswith('H')]
    P_coords = np.array([atom.get_coord() for atom in P_atoms])
    P_coords_pcd = o3d.geometry.PointCloud()
    P_coords_pcd.points = o3d.utility.Vector3dVector(P_coords)
    P_coords_pcd.transform(T)
    for ix, v in enumerate(P_coords_pcd.points):
        P_atoms[ix].set_coord(v)

    if output_path is not None:
        io = PDBIO()
        for atom in Selection.unfold_entities(structure, 'A'):
            if atom.get_name().startswith('H'):
                parent = atom.get_parent()
                parent.detach_child(atom.get_id())
        io.set_structure(structure)
        io.save(output_path)

    return structure
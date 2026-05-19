import json
import os
import numpy as np
import open3d as o3d

from scipy.spatial import cKDTree
from Bio.PDB import PDBParser, PDBIO, Selection, DSSP
from sklearn.cluster import DBSCAN
from packaging import version
from geometry.open3d_import import *
from default_config.masif_opts import masif_opts

if version.parse('0.12.0') <= version.parse(o3d.__version__):
    ICPConvergenceCriteria = o3d.pipelines.registration.ICPConvergenceCriteria
elif version.parse('0.6.0') < version.parse(o3d.__version__):
    ICPConvergenceCriteria = o3d.registration.ICPConvergenceCriteria
else:
    ICPConvergenceCriteria = o3d.ICPConvergenceCriteria

from masif_mimicry.structure.transforms import apply_transform, get_transformed_struct_from_row

def get_features(params: dict, pdb: str, pid: str, source: bool = True, flip_desc: bool = False) -> dict:
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
        mesh = read_triangle_mesh(ply),
        pcd = read_point_cloud(ply),
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
    point_coords: np.ndarray,
    pdb_path: str,
    k: int = 1,
    distance_cutoff: float = 5.0, 
    rsa_cutoff: float = 0.1, 
    exclude_backbone: bool = True,
    dssp_bin_path: str = "/work/lpdi/bin/sequence/dssp"
    ) -> tuple:
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

def get_atom_coords(pdb_path: str, chain: str, residue: int, atom_name: str) -> np.ndarray:
    """Return Nx3 coordinates for all atoms matching chain/residue/atom_name."""
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('target', pdb_path)
    atoms = [
        atom for atom in struct.get_atoms()
        if atom.get_parent().get_id()[1] == residue
        and atom.get_parent().get_parent().get_id() == chain
        and atom.get_id() == atom_name
    ]
    if len(atoms) == 0:
        raise ValueError(
            f'No atom {atom_name} found on chain {chain} residue {residue} in {pdb_path}'
        )
    return np.array([atom.get_coord() for atom in atoms])


def farthest_point_subsample(coords: np.ndarray, num_points: int, seed_idx: int = 0) -> np.ndarray:
    """Farthest-point sampling on discrete 3D points; returns local indices into coords."""
    n = len(coords)
    if n == 0:
        return np.array([], dtype=int)
    k = min(num_points, n)
    selected = [seed_idx]
    if k == 1:
        return np.array(selected, dtype=int)
    min_dists_sq = np.sum((coords - coords[seed_idx]) ** 2, axis=1)
    min_dists_sq[seed_idx] = -1.0
    for _ in range(k - 1):
        next_idx = int(np.argmax(min_dists_sq))
        selected.append(next_idx)
        new_dists_sq = np.sum((coords - coords[next_idx]) ** 2, axis=1)
        min_dists_sq = np.minimum(min_dists_sq, new_dists_sq)
        min_dists_sq[next_idx] = -1.0
    return np.array(selected, dtype=int)


def select_target_sites_by_radius_fps(
    mesh_vertices: np.ndarray,
    pdb_path: str,
    chain: str,
    residue: int,
    atom_name: str,
    radius: float,
    num_points: int,
) -> np.ndarray:
    """
    Select target surface sites within radius (A) of an atom, then FPS to num_points.
    Returns sorted global vertex indices. Warns if fewer than num_points candidates (B1).
    """
    mesh_coords = np.asarray(mesh_vertices)
    atom_coords = get_atom_coords(pdb_path, chain, residue, atom_name)
    dists = np.linalg.norm(mesh_coords[:, None, :] - atom_coords[None, :, :], axis=-1)
    min_dists = dists.min(axis=1)
    candidates = np.where(min_dists <= radius)[0]
    if len(candidates) == 0:
        raise ValueError(
            f'No surface vertices within {radius} A of {atom_name} '
            f'chain {chain} residue {residue} in {pdb_path}'
        )
    if len(candidates) < num_points:
        print(
            f'WARNING: Only {len(candidates)} vertices within {radius} A of '
            f'{atom_name} chain {chain} residue {residue}; using all '
            f'(requested {num_points}).',
            flush=True,
        )
        return np.sort(candidates)
    candidate_coords = mesh_coords[candidates]
    seed_local = int(np.argmin(min_dists[candidates]))
    local_selected = farthest_point_subsample(candidate_coords, num_points, seed_idx=seed_local)
    return np.sort(candidates[local_selected])


TARGET_SITES_MANIFEST = "target_sites.json"


def write_target_vert_files(p2_all_feats, selected_points_idx, target_ppi_id, output_dir, outward_shift=0.25):
    """Write full geodesic patch vertices for each selected target site."""
    vert_dir = os.path.join(output_dir, "target_vert")
    os.makedirs(vert_dir, exist_ok=True)
    pcd_points = np.asarray(p2_all_feats["pcd"].points)
    centers_path = os.path.join(vert_dir, "target.vert")
    with open(centers_path, "w") as out_centers:
        for site_vix in selected_points_idx:
            center = pcd_points[site_vix]
            out_centers.write("{}, {}, {}\n".format(center[0], center[1], center[2]))
            target_patch, _, _ = get_patch_geo(
                p2_all_feats["pcd"],
                p2_all_feats["indices"],
                site_vix,
                p2_all_feats["desc"],
                flip_normals=False,
                outward_shift=outward_shift,
            )
            vert_path = os.path.join(vert_dir, f"{target_ppi_id}_{site_vix}.vert")
            with open(vert_path, "w") as out_patch:
                for point in target_patch.points:
                    out_patch.write("{}, {}, {}\n".format(point[0], point[1], point[2]))
    print(
        f"Wrote {len(selected_points_idx)} patch files and {centers_path} to {vert_dir}",
        flush=True,
    )


def write_target_sites_manifest(target_run_dir, manifest):
    """Write target_sites.json under target_run_dir."""
    path = os.path.join(target_run_dir, TARGET_SITES_MANIFEST)
    with open(path, "w") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return path


def load_target_run_manifest(target_run_dir):
    """
    Load and validate a prepared target run directory.
    Returns (manifest dict, sites as int ndarray).
    """
    target_run_dir = os.path.abspath(os.path.expanduser(target_run_dir))
    manifest_path = os.path.join(target_run_dir, TARGET_SITES_MANIFEST)
    if not os.path.isfile(manifest_path):
        raise FileNotFoundError(
            f"Missing {TARGET_SITES_MANIFEST} in {target_run_dir}. "
            "Run define_target_sites.py first."
        )
    with open(manifest_path) as f:
        manifest = json.load(f)
    vert_dir = os.path.join(target_run_dir, "target_vert")
    if not os.path.isdir(vert_dir):
        raise FileNotFoundError(
            f"Missing target_vert/ in {target_run_dir}. "
            "Run define_target_sites.py first."
        )
    sites = [int(x) for x in manifest["sites"]]
    if len(sites) == 0:
        raise ValueError(f"No target sites listed in {manifest_path}")
    return manifest, np.array(sites, dtype=int)


def write_partner_pdb_for_clashes(params, target_pdb, target_ppi_id, target_run_dir):
    """Write chain-filtered partner PDB once for downstream clash counting."""
    P2_raw_pdb = os.path.join(
        params["masif_target_root"],
        "data_preparation",
        "00-raw_pdbs",
        f"{target_pdb.split('_')[0]}.pdb",
    )
    partner_chain_ids = parse_partner_chain_ids(target_pdb, target_ppi_id)
    partner_suffix = "".join(partner_chain_ids)
    partner_filename = f"{target_pdb.split('_')[0]}_{partner_suffix}.pdb"
    partner_path = os.path.join(target_run_dir, partner_filename)
    filtered = get_filtered_target_structure(P2_raw_pdb, partner_chain_ids, {})
    io = PDBIO()
    io.set_structure(filtered)
    io.save(partner_path)
    return partner_filename, partner_path


def res2surf(point_coords: np.ndarray, pdb_path: str, chain: str, residue: int, atom_name: str, k: int = 1) -> np.ndarray:
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
    atom_coords = get_atom_coords(pdb_path, chain, residue, atom_name)

    # NOTE: this will register the atom/residue to the closest points
    point_cdktree = cKDTree(point_coords)
    d, nearest_points_idx = point_cdktree.query(atom_coords, k=k)
    
    return nearest_points_idx

def select_patches(
    P_all_feats: dict,
    downsample_rate: int = 5, 
    iface_cutoff: float = 0.0, 
    min_patch_num: int = 50, 
    top_iface_percent: float = 0.0, 
    interface_only: bool = False,
    verbose: bool = True
) -> tuple:
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
    
    if top_iface_percent > 0.0:
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

    patch = PointCloud()
    patch.points = Vector3dVector(patch_pts)
    patch.normals = Vector3dVector(patch_nrmls)
    patch_descs = [Feature(), Feature(), Feature()]
    patch_descs[0].data = descriptors[patch_idxs, :].T
    
    return patch, patch_descs, patch_idxs
    
def multidock(
    source_pt, source_pcd, source_patch_idxs, source_descs, 
    target_pt, target_pcd, target_patch_idxs, target_descs, 
    binder_align: bool = False, ransac_skip: bool = False
):
    ransac_radius=1.5
    ransac_iter=10000
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
            # result = registration_ransac_based_on_feature_matching(
            #     source=source_patch, target=target_patch, source_feature=source_patch_descs[0], target_feature=target_patch_descs[0], 
            #     mutual_filter=False, max_correspondence_distance=ransac_radius, 
            #     estimation_method=TransformationEstimationPointToPoint(False), 
            #     ransac_n = 3, checkers=[CorrespondenceCheckerBasedOnEdgeLength(0.9),
            #     CorrespondenceCheckerBasedOnDistance(1.0),
            #     CorrespondenceCheckerBasedOnNormal(np.pi/2)],
            #     criteria=RANSACConvergenceCriteria(max_iteration=ransac_iter, confidence=0.999999),
            #     seed=42
            # )
            result = registration_ransac_based_on_feature_matching(
                source=source_patch, target=target_patch, source_feature=source_patch_descs[0], target_feature=target_patch_descs[0], 
                max_correspondence_distance=ransac_radius,
                estimation_method=TransformationEstimationPointToPoint(False), ransac_n=3,
                checkers=[CorrespondenceCheckerBasedOnEdgeLength(0.9),
                CorrespondenceCheckerBasedOnDistance(1.0),
                CorrespondenceCheckerBasedOnNormal(np.pi/2)],
                criteria=RANSACConvergenceCriteria(ransac_iter, 500)
            )
            init = result.transformation
        else:
            init = np.identity(4)

        result_icp = registration_icp(
            source=source_patch, target=target_patch,
            max_correspondence_distance=1.5, init=init,
            estimation_method=TransformationEstimationPointToPlane(),
            criteria=ICPConvergenceCriteria()
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


def transform_patch_coords(pcd, patch_indices, site, T):
    """Transform geodesic patch vertex coordinates for a surface site."""
    pts = np.asarray(pcd.points)[patch_indices[site]]
    return apply_transform(pts, T)


def compute_descriptor_score(
    P1_patch_coords,
    P1_descs,
    P1_indices,
    P1_site,
    P2_pcd,
    P2_descs,
    P2_indices,
    P2_site,
):
    """
    MaSIF descriptor match score for aligned P1 patch coordinates vs target patch.
    Returns normalized score in [0, 1), or 0.0 when no neighbors within 5 A.
    """
    P2_patch_coords = np.asarray(P2_pcd.points)[P2_indices[P2_site]]

    P2_patch_ckdtree = cKDTree(P2_patch_coords)
    d_nn, r_nn = P2_patch_ckdtree.query(P1_patch_coords)
    neigh = np.where(d_nn <= 5.0)[0]

    if len(neigh) == 0:
        return 0.0

    P1_patch_desc = P1_descs[P1_indices[P1_site]]
    P2_patch_desc = P2_descs[P2_indices[P2_site]]

    desc_dist_score = np.sum(
        np.square(
            (1 / np.sqrt(np.sum(np.square(P1_patch_desc - P2_patch_desc[r_nn]), axis=1)))[neigh]
        )
    )
    return float(np.tanh(desc_dist_score / 80))


def parse_partner_chain_ids(pdb_identifier, target_ppi_id):
    """
    Return partner chain id(s) for clash counting from a MaSIF PDB identifier.

    MaSIF ids use PDB_p1field_p2field where each field is one chain id or several
    concatenated single-character chain ids (e.g. AB = chains A and B as one partner).

    Examples:
      021structure_C_AB, target_ppi_id p1 -> ['A', 'B']  (p2 partner of chain C)
      021structure_C_AB, target_ppi_id p2 -> ['C']       (p1 partner of chains A+B)
    """
    parts = pdb_identifier.split('_')
    if target_ppi_id == 'p1':
        if len(parts) < 3:
            raise ValueError(
                f'Cannot parse partner chains from {pdb_identifier!r} with target_ppi_id p1'
            )
        raw = parts[2]
    else:
        if len(parts) < 2:
            raise ValueError(
                f'Cannot parse partner chains from {pdb_identifier!r} with target_ppi_id p2'
            )
        raw = parts[1]
    return _chain_ids_for_filter(raw)


def _chain_ids_for_filter(target_chain):
    """Normalize chain argument(s) to a list of chain ids for membership tests."""
    if isinstance(target_chain, (list, tuple)):
        return [str(c) for c in target_chain]
    if isinstance(target_chain, str):
        if len(target_chain) == 1:
            return [target_chain]
        # Concatenated ids, e.g. AB -> chains A and B
        return list(target_chain)
    return [str(target_chain)]


def get_filtered_target_structure(target_pdb_path, target_chain, cache):
    """Parse target PDB once and cache chain-filtered structure for clash counting."""
    chain_ids = _chain_ids_for_filter(target_chain)
    key = (os.path.abspath(target_pdb_path), tuple(chain_ids))
    if key not in cache:
        pdb_parser = PDBParser(QUIET=True)
        target_structure = pdb_parser.get_structure('', target_pdb_path)
        chains_to_remove = [
            chain for chain in target_structure.get_chains()
            if chain.id not in chain_ids
        ]
        for chain in chains_to_remove:
            target_structure[0].detach_child(chain.id)
        cache[key] = target_structure
    return cache[key]


def compute_hit_clash_score(
    P1_pdb,
    target_structure,
    descriptor_score,
    ca_clash_threshold=1.0,
    heavy_atom_clash_threshold=5.0,
):
    """
    Count clashes for a transformed source PDB against a cached target structure.
    Returns (ca_clashes, heavy_clashes, final_score).
    final_score is 0.0 when clash thresholds are exceeded, else descriptor_score.
    """
    pdb_parser = PDBParser(QUIET=True)
    source_structure = pdb_parser.get_structure('', P1_pdb)
    ca_clashes, heavy_clashes = count_clashes(source_structure, target_structure, radius=2.0)
    final_score = descriptor_score
    if ca_clashes > ca_clash_threshold or heavy_clashes > heavy_atom_clash_threshold:
        final_score = 0.0
    return ca_clashes, heavy_clashes, final_score


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
    normalized_desc_dist_score = compute_descriptor_score(
        P1_patch_coords, P1_descs, P1_indices, P1_site,
        P2_pcd, P2_descs, P2_indices, P2_site,
    )
    output = [[0, 0], normalized_desc_dist_score]
    
    if compute_clashes:
        assert 'target_structure' in kwargs, 'Please provide target_structure for clash counting.'
        assert 'target_chain' in kwargs, 'Please provide target_chain for clash counting.'
        target_structure = get_filtered_target_structure(
            kwargs['target_structure'], kwargs['target_chain'], kwargs.get('_target_cache', {}),
        )
        ca_clashes, heavy_clashes, final_score = compute_hit_clash_score(
            P1_pdb,
            target_structure,
            normalized_desc_dist_score,
            ca_clash_threshold=kwargs.get('ca_clash_threshold', 1.0),
            heavy_atom_clash_threshold=kwargs.get('heavy_atom_clash_threshold', 5.0),
        )
        output[0][0], output[0][1] = ca_clashes, heavy_clashes
        output[1] = final_score
        return output, target_structure
    else:
        return output, None

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
    if target_ca_coords.size == 0:
        raise ValueError(
            'Target structure has no CA atoms for clash counting; '
            'check partner chain id and raw PDB chain labels.'
        )
    if target_ca_coords.ndim == 1:
        target_ca_coords = target_ca_coords.reshape(-1, 3)
    target_pcd_tree = cKDTree(target_ca_coords)

    d_nn_ca, _ = target_pcd_tree.query(np.asarray(source_ca_coords), k=1, distance_upper_bound=radius)
    clashing_ca = np.sum(d_nn_ca<=radius)
    
    source_atoms = [atom for atom in source_structure.get_atoms() if not atom.get_name().startswith('H')]
    source_coords = np.array([atom.get_coord() for atom in source_atoms])

    d_nn, _ = target_pcd_tree.query(np.asarray(source_coords), k=1, distance_upper_bound=radius)
    clashing = np.sum(d_nn<=radius)
    
    return clashing_ca, clashing

def _backbone_coords(struct, atoms_to_keep=("N", "CA", "C")):
    coords = []
    for atom in struct.get_atoms():
        name = atom.get_name().strip()
        if name in atoms_to_keep:
            coords.append(atom.get_coord())
    if len(coords) == 0:
        raise ValueError("No backbone atoms found in structure")
    return np.stack(coords)


def _pairwise_rmsd(coords1, coords2):
    if coords1.shape != coords2.shape:
        raise ValueError(f"Coordinate shape mismatch: {coords1.shape} vs {coords2.shape}")
    diff = coords1 - coords2
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=-1))))


def structural_clusters(df, rmsd_thresh=5.0, database_dir=None, database_root=None):
    """
    Cluster hits in a mimicry results table by binding mode (RMSD in target frame).

    Expects columns P1_id, flattened_transform.
    Adds cluster_id, cluster_size, cluster_mean_rmsd (appended by caller column order).

    Rows that cannot be transformed or have incompatible backbone size receive
    cluster_id=-1 and NA for cluster_size / cluster_mean_rmsd.
    """
    if database_dir is None:
        database_dir = database_root
    if database_dir is None:
        raise ValueError("database_dir is required for structural_clusters")

    df = df.copy()
    df["cluster_id"] = -1
    df["cluster_size"] = np.nan
    df["cluster_mean_rmsd"] = np.nan

    if len(df) == 0:
        return df

    positions = []
    valid_indices = []
    for idx, row in df.iterrows():
        try:
            struct = get_transformed_struct_from_row(row, database_dir)
            positions.append(_backbone_coords(struct))
            valid_indices.append(idx)
        except Exception:
            continue

    if len(valid_indices) == 0:
        return df

    if len(valid_indices) == 1:
        idx = valid_indices[0]
        df.loc[idx, "cluster_id"] = 0
        df.loc[idx, "cluster_size"] = 1
        df.loc[idx, "cluster_mean_rmsd"] = 0.0
        return df

    ref_shape = positions[0].shape
    filtered_positions = []
    filtered_indices = []
    for pos, idx in zip(positions, valid_indices):
        if pos.shape == ref_shape:
            filtered_positions.append(pos)
            filtered_indices.append(idx)

    if len(filtered_indices) == 0:
        return df

    if len(filtered_indices) == 1:
        idx = filtered_indices[0]
        df.loc[idx, "cluster_id"] = 0
        df.loc[idx, "cluster_size"] = 1
        df.loc[idx, "cluster_mean_rmsd"] = 0.0
        return df

    n = len(filtered_positions)
    rmsd_vals = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            try:
                rmsd_ij = _pairwise_rmsd(filtered_positions[i], filtered_positions[j])
            except ValueError:
                rmsd_ij = np.inf
            rmsd_vals[i, j] = rmsd_ij
            rmsd_vals[j, i] = rmsd_ij

    labels = DBSCAN(eps=rmsd_thresh, min_samples=2, metric="precomputed").fit_predict(rmsd_vals)
    n_outliers = (labels == -1).sum()
    if n_outliers > 0:
        labels[labels == -1] = np.arange(n_outliers) + labels.max() + 1

    current_cluster_label = 0
    for lb in sorted(set(labels)):
        member_mask = labels == lb
        member_pos = np.where(member_mask)[0]
        cluster_size = int(member_mask.sum())
        sub_rmsd = rmsd_vals[member_mask][:, member_mask]
        mean_rmsd_per_member = sub_rmsd.mean(axis=1)
        for local_i, global_i in enumerate(member_pos):
            idx = filtered_indices[global_i]
            df.loc[idx, "cluster_id"] = current_cluster_label
            df.loc[idx, "cluster_size"] = cluster_size
            df.loc[idx, "cluster_mean_rmsd"] = mean_rmsd_per_member[local_i]
        current_cluster_label += 1

    return df


# --- Backward-compatible re-exports (deprecated; use masif_mimicry.*) ---
import warnings


def _deprecated_import(name, module):
    warnings.warn(
        f"utils.{name} is deprecated; import from masif_mimicry.{module}",
        DeprecationWarning,
        stacklevel=2,
    )


def resolve_database_paths(database_dir: str):
    _deprecated_import("resolve_database_paths", "config.paths")
    from masif_mimicry.config.paths import resolve_database_paths as _fn
    return _fn(database_dir)


def set_params(*, database_dir: str, target_preprocess_dir: str, masif_app: str = "ppi_search") -> dict:
    _deprecated_import("set_params", "config.paths")
    from masif_mimicry.config.paths import set_params as _fn
    return _fn(database_dir=database_dir, target_preprocess_dir=target_preprocess_dir, masif_app=masif_app)


def transform_structure(input_path, T, output_path=None):
    _deprecated_import("transform_structure", "structure.transforms")
    from masif_mimicry.structure.transforms import transform_structure as _fn
    return _fn(input_path, T, output_path=output_path)


def seed_pdb_path(p1_id, database_dir):
    _deprecated_import("seed_pdb_path", "structure.transforms")
    from masif_mimicry.structure.transforms import seed_pdb_path as _fn
    return _fn(p1_id, database_dir)


from masif_mimicry.config.paths import (  # noqa: E402
    BENCHMARK_PDB_SUBDIR,
    BENCHMARK_SURF_SUBDIR,
    PRECOMP_12A_SUBDIR,
    PRECOMP_9A_SUBDIR,
)


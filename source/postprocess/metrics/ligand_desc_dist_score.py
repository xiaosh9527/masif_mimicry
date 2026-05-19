from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser

if __name__ == "__main__":
    import sys
    basedir = Path(__file__).resolve().parent.parent.parent
    sys.path.append(str(basedir))

from utils import resolve_database_paths


def pairwise_distances(coord1, coord2):
    return np.linalg.norm(coord1.reshape(-1, 1, 3) - coord2.reshape(1, -1, 3), axis=-1)
    

def desc_dist_score(query_desc, binder_desc, query_surf_coord, binder_surf_coord, query_idx):

    d = pairwise_distances(query_surf_coord[query_idx], binder_surf_coord)
    binder_idx = d.argmin(axis=1)

    desc_dist = np.sqrt(np.sum(np.square(query_desc[query_idx] - binder_desc[binder_idx]), axis=-1))

    # desc dist score: sum over those within 1.5A
    neigh = np.where(d < 1.5)[0]
    desc_dist_score = np.sum(np.square(1 / desc_dist[neigh]))

    return desc_dist_score


def ligand_desc_dist_score(
    target_name,
    binder_name,
    ligand_name,
    ligand_chain,
    binder_transform,
    target_processed_root,
    database_dir,
):
    db_root, db_prep = resolve_database_paths(database_dir)

    # Load surface coordinates
    target_surf_coord = np.stack([
        np.load(Path(target_processed_root, "data_preparation", "04b-precomputation_12A", "precomputation", target_name, f"p1_{dim}.npy"))
        for dim in ['X', 'Y', 'Z']
    ], axis=1)
    binder_surf_coord = np.stack([
        np.load(Path(db_prep, "04b-precomputation_12A", "precomputation", binder_name, f"p1_{dim}.npy"))
        for dim in ['X', 'Y', 'Z']
    ], axis=1)

    # Load descriptors
    target_desc = np.load(Path(target_processed_root, "descriptors", "sc05", "all_feat", target_name, "p1_desc_flipped.npy"))
    binder_desc = np.load(Path(db_root, "descriptors", "sc05", "all_feat", binder_name, "p1_desc_straight.npy"))

    # Transform binder surface point cloud
    R, t = binder_transform[:3, :3], binder_transform[:3, 3]
    binder_surf_coord = binder_surf_coord @ R.T + t[None, :]

    # Extract relevant surface points on target
    target_pdb = Path(target_processed_root, "data_preparation", "01-benchmark_pdbs", f"{target_name}.pdb")  # target preprocess tree
    target_struct = PDBParser(QUIET=True).get_structure("", target_pdb)
    target_atoms = [a for a in target_struct.get_atoms() if a.element != 'H']
    target_coord = np.array([a.get_coord() for a in target_atoms])
    ligand_mask = np.array([(a.parent.resname == ligand_name) and (a.parent.parent.id == ligand_chain) for a in target_atoms])
    dist_to_ligand = pairwise_distances(target_surf_coord, target_coord[ligand_mask]).min(axis=-1)
    dist_to_protein = pairwise_distances(target_surf_coord, target_coord[~ligand_mask]).min(axis=-1)
    ligand_surf_mask = (dist_to_ligand <= dist_to_protein)

    # Find corresponding surface points on binder and calculate the score
    return desc_dist_score(target_desc, binder_desc, target_surf_coord, binder_surf_coord, query_idx=ligand_surf_mask)

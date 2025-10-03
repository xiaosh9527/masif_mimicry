from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser


class RotoTranslation:
    def __init__(self, R, t):
        self.R = R
        self.t = t.reshape(1, 3)

    def __call__(self, coord):
        # return (coord - self.t1) @ self.R.T + self.t2
        return coord @ self.R.T + self.t


def get_transformation(coord, target_coord):
    """
    Compute transformation (R, t) that moves coord to target_coord: target_coord = R @ coord + t.
    Uses the Kabsch algorithm.
    :param coord: numpy array, shape (n, 3)
    :param target_coord: numpy array, shape (n, 3)
    """
    t1 = coord.mean(axis=0)
    t2 = target_coord.mean(axis=0)
    coord = coord - coord.mean(axis=0, keepdims=True)
    target_coord = target_coord - target_coord.mean(axis=0, keepdims=True)

    # https://en.wikipedia.org/wiki/Kabsch_algorithm
    P = target_coord
    Q = coord
    H = P.T @ Q
    U, S, Vh = np.linalg.svd(H)
    d = np.linalg.det(U) * np.linalg.det(Vh)
    R = U @ np.diag([1, 1, d]) @ Vh

    return RotoTranslation(R=R, t=t2 - t1 @ R.T)


def check_transformation(coord, target_coord, transform):
    transformed_coord = transform(coord)
    max_error = np.linalg.norm(transformed_coord - target_coord, axis=-1).max()
    assert max_error < 1e-3, max_error


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


if __name__ == "__main__":
    from argparse import ArgumentParser
    

    parser = ArgumentParser()
    parser.add_argument("--query_processed_root", type=Path)
    parser.add_argument("--query", type=Path)
    parser.add_argument("--ligand", type=str, help="Format: CHAIN_NAME, e.g. B_Y70")
    parser.add_argument("--binder_processed_root", type=Path)
    parser.add_argument("--binder_pdb", type=Path)
    args = parser.parse_args()

    # Load coordinates and descriptors
    binder_name = '_'.join(args.binder_pdb.stem.split("_")[:-1])  # assumes binder is the output of a MaSIF-search
    binder_struct = PDBParser(QUIET=True).get_structure("", args.binder_pdb)
    binder_ca_coord = np.array([res['CA'].get_coord() for res in binder_struct.get_residues()])

    query_surf_coord = np.stack([
        np.load(Path(args.query_processed_root, "data_preparation", "04b-precomputation_12A", "precomputation", args.query, f"p1_{dim}.npy"))
        for dim in ['X', 'Y', 'Z']
    ], axis=1)
    binder_surf_coord = np.stack([
        np.load(Path(args.binder_processed_root, "data_preparation", "04b-precomputation_12A", "precomputation", binder_name, f"p1_{dim}.npy"))
        for dim in ['X', 'Y', 'Z']
    ], axis=1)

    query_desc = np.load(Path(args.query_processed_root, "descriptors", "sc05", "all_feat", args.query, "p1_desc_flipped.npy"))
    binder_desc = np.load(Path(args.binder_processed_root, "descriptors", "sc05", "all_feat", binder_name, "p1_desc_straight.npy"))

    # Rotate surface point cloud if necessary
    binder_ref_struct = PDBParser(QUIET=True).get_structure(
        "", Path(args.binder_processed_root, "data_preparation", "01-benchmark_pdbs", f"{binder_name}.pdb")
    )
    binder_ref_ca_coord = np.array([res['CA'].get_coord() for res in binder_ref_struct.get_residues()])
    transform = get_transformation(binder_ref_ca_coord, binder_ca_coord)
    check_transformation(binder_ref_ca_coord, binder_ca_coord, transform)

    binder_surf_coord = transform(binder_surf_coord)

    # Extract relevant surface points on query
    # Note: Updated the splitting order to CHAIN_NAME (e.g. B_Y70)
    ligand_chain, ligand_name = args.ligand.split("_")
    query_struct = PDBParser(QUIET=True).get_structure(
        "", Path(args.query_processed_root, "data_preparation", "01-benchmark_pdbs", f"{args.query}.pdb")
    )
    query_atoms = [a for a in query_struct.get_atoms() if a.element != 'H']
    query_coord = np.array([a.get_coord() for a in query_atoms])
    ligand_mask = np.array([(a.parent.resname == ligand_name) and (a.parent.parent.id == ligand_chain) for a in query_atoms])
    dist_to_ligand = pairwise_distances(query_surf_coord, query_coord[ligand_mask]).min(axis=-1)
    dist_to_protein = pairwise_distances(query_surf_coord, query_coord[~ligand_mask]).min(axis=-1)
    ligand_surf_mask = (dist_to_ligand <= dist_to_protein)

    # Find corresponding surface points on binder and calculate the score
    print(desc_dist_score(query_desc, binder_desc, query_surf_coord, binder_surf_coord, query_idx=ligand_surf_mask))
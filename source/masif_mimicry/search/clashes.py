"""Bio.PDB clash counting for search-time hit scoring (not postprocess RDKit clashes)."""

import numpy as np
from scipy.spatial import cKDTree


def count_clashes(source_structure, target_structure, radius=2.0):
    """
    Count clashing CA and heavy atoms between Bio.PDB structures (search path).
    """
    source_ca_coords = np.array(
        [atom.get_coord() for atom in source_structure.get_atoms() if atom.get_id() == "CA"]
    )

    target_ca_coords = np.array(
        [atom.get_coord() for atom in target_structure.get_atoms() if atom.get_id() == "CA"]
    )
    if target_ca_coords.size == 0:
        raise ValueError(
            "Target structure has no CA atoms for clash counting; "
            "check partner chain id and raw PDB chain labels."
        )
    if target_ca_coords.ndim == 1:
        target_ca_coords = target_ca_coords.reshape(-1, 3)
    target_pcd_tree = cKDTree(target_ca_coords)

    d_nn_ca, _ = target_pcd_tree.query(
        np.asarray(source_ca_coords), k=1, distance_upper_bound=radius
    )
    clashing_ca = np.sum(d_nn_ca <= radius)

    source_atoms = [
        atom for atom in source_structure.get_atoms() if not atom.get_name().startswith("H")
    ]
    source_coords = np.array([atom.get_coord() for atom in source_atoms])

    d_nn, _ = target_pcd_tree.query(
        np.asarray(source_coords), k=1, distance_upper_bound=radius
    )
    clashing = np.sum(d_nn <= radius)

    return clashing_ca, clashing

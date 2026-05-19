import os

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Selection
from geometry.open3d_import import PointCloud, Vector3dVector

from masif_mimicry.config.paths import BENCHMARK_PDB_SUBDIR, resolve_database_paths


def apply_transform(coords, T):
    """Apply a 4x4 homogeneous transform to Nx3 coordinates."""
    coords = np.asarray(coords, dtype=float)
    if coords.ndim == 1:
        coords = coords.reshape(1, -1)
    n = coords.shape[0]
    coords_h = np.hstack([coords, np.ones((n, 1))])
    return (np.asarray(T) @ coords_h.T).T[:, :3]


def transform_structure(input_path, T, output_path=None):
    """
    Transform a PDB structure by applying transformation matrix T (Open3D convention).

    Returns the transformed Bio.PDB Structure. Optionally writes output_path.
    """
    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure("", input_path)
    p_atoms = [atom for atom in structure.get_atoms() if not atom.get_name().startswith("H")]
    p_coords = np.array([atom.get_coord() for atom in p_atoms])
    p_coords_pcd = PointCloud()
    p_coords_pcd.points = Vector3dVector(p_coords)
    p_coords_pcd.transform(T)
    for ix, v in enumerate(p_coords_pcd.points):
        p_atoms[ix].set_coord(v)

    if output_path is not None:
        io = PDBIO()
        for atom in Selection.unfold_entities(structure, "A"):
            if atom.get_name().startswith("H"):
                parent = atom.get_parent()
                parent.detach_child(atom.get_id())
        io.set_structure(structure)
        io.save(output_path)

    return structure


def seed_pdb_path(p1_id, database_dir):
    """Resolve preprocessed seed PDB under database_dir/data_preparation/01-benchmark_pdbs/."""
    _, db_prep = resolve_database_paths(database_dir)
    return os.path.join(db_prep, BENCHMARK_PDB_SUBDIR, f"{p1_id}.pdb")


def get_transformed_struct_from_row(row, database_dir):
    """Apply flattened_transform to the seed PDB using mimicry/Open3D conventions."""
    pdb_path = seed_pdb_path(row.P1_id, database_dir)
    t = np.array(list(map(float, row.flattened_transform.split(",")))).reshape(4, 4)
    return transform_structure(pdb_path, t, output_path=None)

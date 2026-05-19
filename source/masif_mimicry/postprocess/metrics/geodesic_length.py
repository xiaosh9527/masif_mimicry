"""
Code from Yanxiang Meng
"""

import argparse
import numpy as np
import networkx as nx
from skimage.morphology import skeletonize
from scipy.ndimage import median_filter
from Bio.PDB import PDBParser

"""
Absolute length: Longest path length in the skeletonized voxel grid.
Normalized length: Absolute length divided by the cube root of the object volume.
The calculations are based on the Van-der-Waals radii of atoms, which are used to
create a voxel grid representation of the molecular structure.
"""

# Van-der-Waals radii (Å)
VDW_RADII = {
    "H": 1.2, "HE": 1.4, "LI": 1.82, "BE": 1.53, "B": 1.92, "C": 1.70, "N": 1.55,
    "O": 1.52, "F": 1.47, "NE": 1.54, "NA": 2.27, "MG": 1.73, "AL": 1.84, "SI": 2.10,
    "P": 1.80, "S": 1.80, "CL": 1.75, "AR": 1.88,
}


def load_structure(input_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("mol", input_path)
    atoms = []
    for atom in structure.get_atoms():
        element = atom.element.strip().upper() if atom.element else atom.get_name()[0]
        radius = VDW_RADII.get(element, 1.7)
        atoms.append((atom.get_coord(), radius))
    return atoms


def create_voxel_grid(atoms, *, spacing: float = 1.0, margin: float = 5.0):
    positions = np.array([pos for pos, _ in atoms])
    min_c = positions.min(axis=0) - margin
    max_c = positions.max(axis=0) + margin

    dims = np.ceil((max_c - min_c) / spacing).astype(int) + 1
    grid = np.zeros(dims, dtype=bool)

    for pos, radius in atoms:
        centre = ((pos - min_c) / spacing).astype(int)
        r_vox = int(np.ceil(radius / spacing))
        xmin, xmax = max(centre[0] - r_vox, 0), min(centre[0] + r_vox + 1, dims[0])
        ymin, ymax = max(centre[1] - r_vox, 0), min(centre[1] + r_vox + 1, dims[1])
        zmin, zmax = max(centre[2] - r_vox, 0), min(centre[2] + r_vox + 1, dims[2])
        for i in range(xmin, xmax):
            for j in range(ymin, ymax):
                for k in range(zmin, zmax):
                    v_c = min_c + np.array([i, j, k]) * spacing
                    if np.sum((v_c - pos) ** 2) <= radius ** 2:
                        grid[i, j, k] = True
    return grid, min_c


def longest_path_skeleton(skel):
    voxels = np.argwhere(skel)
    G = nx.Graph()
    for v in map(tuple, voxels):
        G.add_node(v)
    for v in map(tuple, voxels):
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                for dz in [-1, 0, 1]:
                    if dx == dy == dz == 0:
                        continue
                    n = (v[0] + dx, v[1] + dy, v[2] + dz)
                    if n in G:
                        G.add_edge(v, n, weight=np.linalg.norm(np.array(v) - np.array(n)))
    endpoints = [n for n, deg in G.degree() if deg == 1] or list(G.nodes)
    source = endpoints[0]
    dists = nx.single_source_dijkstra_path_length(G, source, weight="weight")
    farthest = max(dists, key=dists.get)
    dists2, paths = nx.single_source_dijkstra(G, farthest, weight="weight")
    other = max(dists2, key=dists2.get)
    return paths[other], dists2[other]


def compute_geodesic_lengths(input_path, *, return_abs: bool = True, return_norm: bool = True, print_results: bool = False):
    atoms = load_structure(str(input_path))
    spacing = 1.0
    voxel_grid, _ = create_voxel_grid(atoms, spacing=spacing)

    smoothed = median_filter(voxel_grid.astype(float), size=3) > 0.5
    _, abs_len = longest_path_skeleton(skeletonize(smoothed))
    volume = np.sum(voxel_grid) * (spacing ** 3)
    norm_len = abs_len / (volume ** (1 / 3)) if abs_len > 0 else 0

    if print_results:
        if return_abs:
            print(abs_len)
        if return_norm:
            print(norm_len)

    return (abs_len if return_abs else None, norm_len if return_norm else None)


if __name__ == "__main__":
    argp = argparse.ArgumentParser(description="Geodesic length calculation")
    argp.add_argument("--input", required=True)
    argp.add_argument("--abs", dest="abs_flag", action="store_true")
    argp.add_argument("--norm", dest="norm_flag", action="store_true")
    cli_args = argp.parse_args()

    if not cli_args.abs_flag and not cli_args.norm_flag:
        cli_args.abs_flag = cli_args.norm_flag = True

    abs_l, norm_l = compute_geodesic_lengths(
        cli_args.input,
        return_abs=cli_args.abs_flag,
        return_norm=cli_args.norm_flag,
        print_results=True,
    )
#!/usr/bin/env python3
import argparse
import numpy as np
import networkx as nx
from skimage.morphology import skeletonize_3d
from scipy.ndimage import median_filter
from Bio.PDB import PDBParser

"""
geodesic_length.py
Yanxiang Meng - 14.03.2025

Calculates geosedic length to meausre how elongated an object is.
norm_geosedic_length = abs_geodesic_length / (object_volume ** (1/3)); effectively the ratio between the longest dimension / average dimension
"""

# Van der Waals radii (in Å) for some common elements.
VDW_RADII = {
    "H": 1.2,  "HE": 1.4, "LI": 1.82, "BE": 1.53,
    "B": 1.92, "C": 1.70, "N": 1.55,  "O": 1.52,
    "F": 1.47, "NE": 1.54, "NA": 2.27, "MG": 1.73,
    "AL": 1.84, "SI": 2.10, "P": 1.80,  "S": 1.80,
    "CL": 1.75, "AR": 1.88,
    # Extend as needed...
}

def load_structure(input_path):
    """
    Load a PDB structure using BioPython.
    Returns a list of tuples (position, radius) for each atom.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("mol", input_path)
    
    atoms_data = []
    for atom in structure.get_atoms():
        element = atom.element.strip().upper() if atom.element is not None else atom.get_name()[0]
        radius = VDW_RADII.get(element, 1.7)  # default radius if element not found
        pos = atom.get_coord()
        atoms_data.append((pos, radius))
    return atoms_data

def create_voxel_grid(atoms, spacing=1.0, margin=5.0):
    """
    Create a 3D binary voxel grid from a list of atoms.
    Each atom is represented as a sphere with its own radius.
    """
    positions = np.array([atom[0] for atom in atoms])
    min_coords = positions.min(axis=0) - margin
    max_coords = positions.max(axis=0) + margin
    grid_dims = np.ceil((max_coords - min_coords) / spacing).astype(int) + 1
    grid = np.zeros(grid_dims, dtype=bool)
    
    for pos, radius in atoms:
        center = ((pos - min_coords) / spacing).astype(int)
        r_vox = int(np.ceil(radius / spacing))
        xmin = max(center[0] - r_vox, 0)
        xmax = min(center[0] + r_vox + 1, grid_dims[0])
        ymin = max(center[1] - r_vox, 0)
        ymax = min(center[1] + r_vox + 1, grid_dims[1])
        zmin = max(center[2] - r_vox, 0)
        zmax = min(center[2] + r_vox + 1, grid_dims[2])
        for i in range(xmin, xmax):
            for j in range(ymin, ymax):
                for k in range(zmin, zmax):
                    voxel_center = min_coords + np.array([i, j, k]) * spacing
                    if np.sum((voxel_center - pos)**2) <= radius**2:
                        grid[i, j, k] = True
    return grid, min_coords

def longest_path_skeleton(skel):
    """
    Build a graph from the skeleton (binary 3D array) and compute its longest geodesic path.
    Returns a list of voxel indices (as tuples) along the longest path and its total length.
    """
    skel_voxels = np.argwhere(skel)
    G = nx.Graph()
    for voxel in map(tuple, skel_voxels):
        G.add_node(voxel)
    for voxel in map(tuple, skel_voxels):
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                for dz in [-1, 0, 1]:
                    if dx == 0 and dy == 0 and dz == 0:
                        continue
                    neighbor = (voxel[0] + dx, voxel[1] + dy, voxel[2] + dz)
                    if neighbor in G:
                        weight = np.linalg.norm(np.array(voxel) - np.array(neighbor))
                        G.add_edge(voxel, neighbor, weight=weight)
    endpoints = [node for node, degree in G.degree() if degree == 1]
    if not endpoints:
        endpoints = list(G.nodes)
    source = endpoints[0]
    distances = nx.single_source_dijkstra_path_length(G, source, weight='weight')
    farthest_node = max(distances, key=distances.get)
    distances2, paths = nx.single_source_dijkstra(G, farthest_node, weight='weight')
    other_end = max(distances2, key=distances2.get)
    longest_path = paths[other_end]
    return longest_path, distances2[other_end]

def resample_path(world_coords, interval=1.0):
    """
    Resample a list of 3D points (world coordinates) so that points are approximately 'interval' apart.
    """
    if not world_coords:
        return []
    sampled = [world_coords[0]]
    accumulated = 0.0
    prev = world_coords[0]
    for pt in world_coords[1:]:
        d = np.linalg.norm(pt - prev)
        accumulated += d
        if accumulated >= interval:
            sampled.append(pt)
            accumulated = 0.0
        prev = pt
    return sampled

def main():
    parser = argparse.ArgumentParser(
        description="Compute geodesic length of a PDB model via skeletonization of its 3D object."
    )
    parser.add_argument('--input', required=True, help="Path to input PDB file")
    parser.add_argument('--output', required=False, help="Path to output file for sampled skeleton vertices")
    parser.add_argument('--abs', dest='abs_flag', action='store_true',
                        help="Return the absolute geodesic length as a numeric value")
    parser.add_argument('--norm', dest='norm_flag', action='store_true',
                        help="Return the normalized geodesic length (absolute length / volume) as a numeric value")
    args = parser.parse_args()
    
    atoms = load_structure(args.input)
    spacing = 1.0  # voxel size in Å
    voxel_grid, origin = create_voxel_grid(atoms, spacing=spacing, margin=5.0)
    
    # Apply median filter for smoothing the voxel grid.
    # Since the grid is boolean, convert to float (0.0, 1.0) for filtering.
    float_grid = voxel_grid.astype(float)
    size = 3  # Adjust size as needed for smoothing
    smoothed = median_filter(float_grid, size=size)
    # Threshold back to binary
    smoothed_binary = smoothed > 0.5
    
    object_volume = np.sum(voxel_grid) * (spacing**3)
    
    skeleton = skeletonize_3d(smoothed_binary)
    longest_voxel_path, geodesic_length = longest_path_skeleton(skeleton)
    normalized_length = geodesic_length / (object_volume ** (1/3)) if geodesic_length > 0 else 0
    
    # Output results based on provided flags.
    if args.abs_flag and not args.norm_flag:
        print(geodesic_length)
    elif args.norm_flag and not args.abs_flag:
        print(normalized_length)
    elif args.abs_flag and args.norm_flag:
        print(geodesic_length)
        print(normalized_length)
    else:
        print("Absolute geodesic length:", geodesic_length)
        print("Normalized geodesic length:", normalized_length)
    
    if args.output:
        world_path = [np.array(idx) * spacing + origin for idx in longest_voxel_path]
        sampled_points = resample_path(world_path, interval=1.0)
        with open(args.output, 'w') as f:
            for pt in sampled_points:
                f.write("{},{},{}\n".format(pt[0], pt[1], pt[2]))

if __name__ == "__main__":
    main()
import json
import os

import numpy as np
from Bio.PDB import DSSP, PDBIO, PDBParser
from scipy.spatial import cKDTree

from masif_mimicry.search.features import get_patch_geo
from masif_mimicry.search.manifest import TARGET_SITES_MANIFEST
from masif_mimicry.utils.parse_structures import (
    get_filtered_target_structure,
    parse_partner_chain_ids,
)


def surf2atom(
    point_coords: np.ndarray,
    pdb_path: str,
    k: int = 1,
    distance_cutoff: float = 5.0,
    rsa_cutoff: float = 0.1,
    exclude_backbone: bool = True,
    dssp_bin_path: str = "/work/lpdi/bin/sequence/dssp",
) -> tuple:
    """Map surface points to closest surface-exposed atoms via DSSP."""
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("", pdb_path)

    dssp = DSSP(struct[0], pdb_path, dssp=dssp_bin_path)
    ss_atom_list = np.concatenate(
        [[x[2]] * len(list(y.get_atoms())) for x, y in zip(dssp, struct.get_residues())]
    )
    rsa_atom_list = np.concatenate(
        [[x[3]] * len(list(y.get_atoms())) for x, y in zip(dssp, struct.get_residues())]
    )
    rsa_atom_list = np.array([float(x) if x != "NA" else 1.0 for x in rsa_atom_list])

    if rsa_cutoff is not None:
        if exclude_backbone:
            atoms = np.array(
                [
                    a
                    for a, b in zip(struct.get_atoms(), rsa_atom_list)
                    if b >= rsa_cutoff and a.get_id() not in ["N", "CA", "C", "O"]
                ]
            )
            ss_label = np.array(
                [
                    c
                    for a, b, c in zip(struct.get_atoms(), rsa_atom_list, ss_atom_list)
                    if b >= rsa_cutoff and a.get_id() not in ["N", "CA", "C", "O"]
                ]
            )
        else:
            atoms = np.array([a for a, b in zip(struct.get_atoms(), rsa_atom_list) if b >= rsa_cutoff])
            ss_label = np.array([a for a, b in zip(ss_atom_list, rsa_atom_list) if b >= rsa_cutoff])
    else:
        if exclude_backbone:
            atoms = np.array([a for a in struct.get_atoms() if a.get_id() not in ["N", "CA", "C", "O"]])
            ss_label = np.array(
                [a for a, b in zip(ss_atom_list, struct.get_atoms()) if b.get_id() not in ["N", "CA", "C", "O"]]
            )
        else:
            atoms = np.array(list(struct.get_atoms()))
            ss_label = ss_atom_list

    atom_coords = np.array([atom.get_coord() for atom in atoms])
    atom_cdktree = cKDTree(atom_coords)
    d, nearest_atoms_idx = atom_cdktree.query(point_coords, k=k)

    if k == 1:
        nearest_atom = [a for a in atoms[nearest_atoms_idx]]
        nearest_atom_ss = [a for a in ss_label[nearest_atoms_idx]]
    else:
        nearest_atom = []
        nearest_atom_ss = []
        for i in range(d.shape[0]):
            nearest_atom.append(
                [atoms[nearest_atoms_idx[i, j]] for j in range(d.shape[1]) if d[i, j] <= distance_cutoff]
            )
            nearest_atom_ss.append(
                [ss_label[nearest_atoms_idx[i, j]] for j in range(d.shape[1]) if d[i, j] <= distance_cutoff]
            )

    return nearest_atom, nearest_atom_ss


def get_atom_coords(pdb_path: str, chain: str, residue: int, atom_name: str) -> np.ndarray:
    """Return Nx3 coordinates for all atoms matching chain/residue/atom_name."""
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("target", pdb_path)
    atoms = [
        atom
        for atom in struct.get_atoms()
        if atom.get_parent().get_id()[1] == residue
        and atom.get_parent().get_parent().get_id() == chain
        and atom.get_id() == atom_name
    ]
    if len(atoms) == 0:
        raise ValueError(
            f"No atom {atom_name} found on chain {chain} residue {residue} in {pdb_path}"
        )
    return np.array([atom.get_coord() for atom in atoms])


def farthest_point_subsample(coords: np.ndarray, num_points: int, seed_idx: int = 0) -> np.ndarray:
    """Farthest-point sampling; returns local indices into coords."""
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
    """Select target sites within radius of an atom, then FPS subsample."""
    mesh_coords = np.asarray(mesh_vertices)
    atom_coords = get_atom_coords(pdb_path, chain, residue, atom_name)
    dists = np.linalg.norm(mesh_coords[:, None, :] - atom_coords[None, :, :], axis=-1)
    min_dists = dists.min(axis=1)
    candidates = np.where(min_dists <= radius)[0]
    if len(candidates) == 0:
        raise ValueError(
            f"No surface vertices within {radius} A of {atom_name} "
            f"chain {chain} residue {residue} in {pdb_path}"
        )
    if len(candidates) < num_points:
        print(
            f"WARNING: Only {len(candidates)} vertices within {radius} A of "
            f"{atom_name} chain {chain} residue {residue}; using all "
            f"(requested {num_points}).",
            flush=True,
        )
        return np.sort(candidates)
    candidate_coords = mesh_coords[candidates]
    seed_local = int(np.argmin(min_dists[candidates]))
    local_selected = farthest_point_subsample(candidate_coords, num_points, seed_idx=seed_local)
    return np.sort(candidates[local_selected])


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
    """Load and validate a prepared target run directory."""
    target_run_dir = os.path.abspath(os.path.expanduser(target_run_dir))
    manifest_path = os.path.join(target_run_dir, TARGET_SITES_MANIFEST)
    if not os.path.isfile(manifest_path):
        raise FileNotFoundError(
            f"Missing {TARGET_SITES_MANIFEST} in {target_run_dir}. "
            "Run python -m masif_mimicry.define_target_sites first."
        )
    with open(manifest_path) as f:
        manifest = json.load(f)
    vert_dir = os.path.join(target_run_dir, "target_vert")
    if not os.path.isdir(vert_dir):
        raise FileNotFoundError(
            f"Missing target_vert/ in {target_run_dir}. "
            "Run python -m masif_mimicry.define_target_sites first."
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
    """Return indices of closest surface points to a residue/atom."""
    atom_coords = get_atom_coords(pdb_path, chain, residue, atom_name)
    point_cdktree = cKDTree(point_coords)
    d, nearest_points_idx = point_cdktree.query(atom_coords, k=k)
    return nearest_points_idx

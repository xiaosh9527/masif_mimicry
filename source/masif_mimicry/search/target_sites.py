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


def parse_site_idx_list(site_idx: str) -> list:
    """Parse comma-separated MaSIF surface site indices (e.g. '1207,940,300')."""
    if not site_idx or not str(site_idx).strip():
        raise ValueError("--site_idx must be a non-empty comma-separated list of integers")
    out = []
    seen = set()
    for part in str(site_idx).split(","):
        part = part.strip()
        if not part:
            continue
        try:
            idx = int(part)
        except ValueError as exc:
            raise ValueError(f"Invalid site index {part!r} in --site_idx (expected integers)") from exc
        if idx in seen:
            continue
        seen.add(idx)
        out.append(idx)
    if not out:
        raise ValueError("--site_idx contained no valid integers")
    return out


def validate_site_indices(site_indices, num_sites: int, pdb_identifier: str) -> np.ndarray:
    """Ensure manual site indices are in range for this target's descriptor table."""
    arr = np.asarray(site_indices, dtype=int)
    bad = arr[(arr < 0) | (arr >= num_sites)]
    if len(bad) > 0:
        raise ValueError(
            f"Site index(es) {bad.tolist()} out of range [0, {num_sites - 1}] for {pdb_identifier} "
            f"({num_sites} MaSIF patch centers)"
        )
    return arr


def poisson_disk_subsample(mesh, num_points: int) -> np.ndarray:
    """Poisson-disk sample mesh vertices; returns sorted vertex indices (non-deterministic)."""
    if not mesh.has_vertex_normals():
        mesh.compute_vertex_normals()
    sampled_points = mesh.sample_points_poisson_disk(num_points)
    mesh_coords = np.asarray(mesh.vertices)
    squared_dists = np.sum(
        np.square(
            np.asarray(sampled_points.points).reshape(-1, 1, 3)
            - mesh_coords.reshape(1, -1, 3)
        ),
        axis=-1,
    )
    subsampled_indices = np.argmin(squared_dists, axis=-1)
    return np.sort(subsampled_indices)


def _resolve_residue_id(structure, chain: str, residue_number: int):
    """Return the BioPython residue id tuple for a unique residue number on chain."""
    matching = [
        res.id
        for res in structure[0][chain].get_residues()
        if res.id[1] == residue_number
    ]
    if len(matching) == 0:
        raise ValueError(
            f"No residue {residue_number} on chain {chain}"
        )
    if len(matching) > 1:
        raise ValueError(
            f"Residue number {residue_number} on chain {chain} is not unique: {matching}"
        )
    return matching[0]


def select_target_sites_by_grid(
    p2_all_feats: dict,
    query_pdb: str,
    chain: str,
    residue: int,
    num_points: int,
    distance_cutoff: float = 4.0,
) -> np.ndarray:
    """
    Grid target site selection (tricomplex search_grid.py logic).

    Locate a residue in an arbitrary query PDB, find MaSIF mesh vertices within
    distance_cutoff of any heavy atom, Poisson-disk subsample the mesh, and keep
    the num_points vertices closest to the residue.
    """
    mesh = p2_all_feats["mesh"]
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("query", query_pdb)
    if chain not in structure[0]:
        raise ValueError(f"Chain {chain!r} not found in {query_pdb}")
    resid = _resolve_residue_id(structure, chain, residue)

    residue_coords = np.stack(
        [
            atom.get_coord()
            for atom in structure[0][chain][resid].get_atoms()
            if atom.element != "H"
        ]
    )
    if len(residue_coords) == 0:
        raise ValueError(
            f"No heavy atoms for chain {chain} residue {residue} in {query_pdb}"
        )

    mesh_coords = np.asarray(mesh.vertices)
    dists = np.sqrt(
        np.sum(
            np.square(mesh_coords.reshape(-1, 1, 3) - residue_coords.reshape(1, -1, 3)),
            axis=-1,
        )
    )
    residue_surface_indices = np.where(np.any(dists < distance_cutoff, axis=-1))[0]
    if len(residue_surface_indices) == 0:
        raise ValueError(
            f"No mesh vertices within {distance_cutoff} A of chain {chain} "
            f"residue {residue} in {query_pdb}"
        )

    num_sites = min(num_points, len(residue_surface_indices))
    poisson_n = int(num_sites * len(mesh_coords) / len(residue_surface_indices))
    poisson_n = max(poisson_n, num_sites)
    mesh_inds = poisson_disk_subsample(mesh, poisson_n)
    closest_to_residue = np.argsort(np.min(dists[mesh_inds, :], axis=-1))
    target_vertices = mesh_inds[closest_to_residue[:num_sites]]
    return np.sort(target_vertices)


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
    partner_chain_ids = parse_partner_chain_ids(target_pdb, target_ppi_id)
    if partner_chain_ids is None:
        return None, None

    P2_raw_pdb = os.path.join(
        params["masif_target_root"],
        "data_preparation",
        "00-raw_pdbs",
        f"{target_pdb.split('_')[0]}.pdb",
    )
    partner_suffix = "".join(partner_chain_ids)
    partner_filename = f"{target_pdb.split('_')[0]}_{partner_suffix}.pdb"
    partner_path = os.path.join(target_run_dir, partner_filename)
    filtered = get_filtered_target_structure(P2_raw_pdb, partner_chain_ids, {})
    io = PDBIO()
    io.set_structure(filtered)
    io.save(partner_path)
    return partner_filename, partner_path



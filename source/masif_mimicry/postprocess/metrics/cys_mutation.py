import numpy as np

from masif_mimicry.postprocess.structures import maybe_load_structure
from masif_mimicry.postprocess.metrics.sasa import extract_ligand_residue


def find_closest_c_beta(query, target_pdb):

    query_struct = maybe_load_structure(query["pdb"])[0]
    ligand = extract_ligand_residue(query_struct, query["chain"], query["resname"])
    query_coord = ligand[query["atom_name"]].get_coord()[None, :]

    target_struct = maybe_load_structure(target_pdb)[0]
    c_betas = [atom for atom in target_struct.get_atoms() if atom.name == "CB"]

    c_beta_coords = np.stack([a.get_coord() for a in c_betas])
    dists = np.sqrt(np.sum((query_coord - c_beta_coords)**2, axis=1))
    min_idx = np.argmin(dists)

    closest_c_beta = c_betas[min_idx]
    closest_c_beta_id = f"{closest_c_beta.parent.parent.id}_{closest_c_beta.parent.get_resname()}{closest_c_beta.parent.id[1]}"

    return {
        "closest_c_beta": closest_c_beta_id,
        "closest_c_beta_distance": dists[min_idx],
    }


def find_smallest_possible_disulfide_length_after_cys_mutation(
    sulfur_coord,
    target_residue,
):
    """
    Assuming `target_residue` is mutated to a cysteine, compute the smallest
    possible S-S distance if a disulfide bond is formed with a sulfur atom
    located in `sulfur_coord`.
    """
    # CB_S_BOND_LENGTH = 1.83  # Angstrom
    # CA_CB_S_ANGLE = ?

    # From PDB 8G9Q
    CB_S_BOND_LENGTH = 1.8  # Angstrom
    CA_CB_S_ANGLE = 117.1 / 180 * np.pi

    bond_axis = target_residue['CB'].get_coord() - target_residue['CA'].get_coord()
    bond_axis = bond_axis / np.linalg.norm(bond_axis)  # normal vector of the plane the sulfur atom sits on

    r = abs(np.sin(np.pi - CA_CB_S_ANGLE) * CB_S_BOND_LENGTH)
    center_coord = target_residue['CB'].get_coord() + np.cos(np.pi - CA_CB_S_ANGLE) * CB_S_BOND_LENGTH * bond_axis

    delta = sulfur_coord - center_coord
    sulfur_proj = delta - (delta @ bond_axis) * bond_axis
    sulfur_proj = center_coord + sulfur_proj

    d_prime = abs(np.linalg.norm(sulfur_proj - center_coord) - r)
    h = abs(delta @ bond_axis)

    # # Brute force for debugging
    # x_vec = (sulfur_proj - center_coord) / np.linalg.norm((sulfur_proj - center_coord) )
    # y_vec = np.cross(bond_axis, x_vec) / np.linalg.norm(np.cross(bond_axis, x_vec))
    # dists = []
    # for phi in np.linspace(0, 2 * np.pi, 360 * 2):
    #     pos = r * np.sin(phi) * x_vec + r * np.cos(phi) * y_vec + center_coord
    #     dists.append(np.linalg.norm(sulfur_coord - pos))
    # print(min(dists), max(dists))

    return np.sqrt(h**2 + d_prime**2)
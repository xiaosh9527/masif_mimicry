import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Polypeptide import is_aa


def find_iface(target_struct: Structure, binder_struct: Structure, max_heavy_atom_dist: float = 5.0, target_patch_coord: np.array = None):
    target_atoms = [a for a in target_struct.get_atoms() if a.element != 'H' and is_aa(a.parent, standard=True)]
    binder_atoms = [a for a in binder_struct.get_atoms() if a.element != 'H' and is_aa(a.parent, standard=True)]

    if target_patch_coord is not None:
        # Keep only atoms close to the selected target patch
        atom_filter = lambda atom: np.linalg.norm(atom.get_coord()[None, :] - target_patch_coord, axis=-1).min() <= 4.5
        target_atoms = [a for a in target_atoms if atom_filter(a)]

    target_resi = np.array([a.parent.id[1] for a in target_atoms])
    binder_resi = np.array([a.parent.id[1] for a in binder_atoms])

    target_coord = np.stack([a.get_coord() for a in target_atoms])
    binder_coord = np.stack([a.get_coord() for a in binder_atoms])

    dists = np.linalg.norm(target_coord[:, None, :] - binder_coord[None, :, :], axis=-1)
    interacting = (dists <= max_heavy_atom_dist)

    target_iface_resi = np.unique(target_resi[interacting.any(axis=1)]).tolist()
    binder_iface_resi = np.unique(binder_resi[interacting.any(axis=0)]).tolist()

    return target_iface_resi, binder_iface_resi


def compute_avg_bfactor(binder_struct, iface_residues):
    if len(iface_residues) < 1:
        return None
    iface_residues = [res for res in binder_struct.get_residues() if res.id[1] in iface_residues]
    bfactors = [res['CA'].get_bfactor() for res in iface_residues]
    return np.mean(bfactors)


def compute_binder_interface_metrics(
        binder_struct,
        target_struct,
        target_patch_coord=None,
    ):
    """
    Compute binder interface metrics including:
      - Number of interface residues and average CA B-factor.
      - Interface residue numbering.
    Returns a dictionary with the computed metrics.
    """
    metrics = {}

    try:
        target_iface, matched_iface = find_iface(target_struct, binder_struct, target_patch_coord=target_patch_coord)
        metrics['target_iface_resi'] = ','.join(map(str, target_iface))
        metrics['matched_iface_resi'] = ','.join(map(str, matched_iface))

        metrics['matched_iface_n_resi'] = len(matched_iface)
        metrics['matched_iface_plddt'] = compute_avg_bfactor(binder_struct, matched_iface)
    except Exception as e:
        print(f"Error computing binder interface metrics: {e}")
        metrics['target_iface_resi'] = None
        metrics['matched_iface_resi'] = None
        metrics['matched_iface_n_resi'] = None
        metrics['matched_iface_plddt'] = None
        return metrics

    return metrics

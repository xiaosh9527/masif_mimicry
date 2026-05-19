import numpy as np
from Bio.PDB import PDBParser
from scipy.spatial import cKDTree

from masif_mimicry.search.clashes import count_clashes
from masif_mimicry.utils.parse_structures import get_filtered_target_structure


def compute_descriptor_score(
    P1_patch_coords,
    P1_descs,
    P1_indices,
    P1_site,
    P2_pcd,
    P2_descs,
    P2_indices,
    P2_site,
):
    """MaSIF descriptor match score for aligned P1 patch vs target patch."""
    P2_patch_coords = np.asarray(P2_pcd.points)[P2_indices[P2_site]]

    P2_patch_ckdtree = cKDTree(P2_patch_coords)
    d_nn, r_nn = P2_patch_ckdtree.query(P1_patch_coords)
    neigh = np.where(d_nn <= 5.0)[0]

    if len(neigh) == 0:
        return 0.0

    P1_patch_desc = P1_descs[P1_indices[P1_site]]
    P2_patch_desc = P2_descs[P2_indices[P2_site]]

    desc_dist_score = np.sum(
        np.square(
            (1 / np.sqrt(np.sum(np.square(P1_patch_desc - P2_patch_desc[r_nn]), axis=1)))[neigh]
        )
    )
    return float(np.tanh(desc_dist_score / 80))


def compute_hit_clash_score(
    P1_pdb,
    target_structure,
    descriptor_score,
    ca_clash_threshold=1.0,
    heavy_atom_clash_threshold=5.0,
):
    """Count clashes for a transformed source PDB; return score zeroed if over threshold."""
    pdb_parser = PDBParser(QUIET=True)
    source_structure = pdb_parser.get_structure("", P1_pdb)
    ca_clashes, heavy_clashes = count_clashes(source_structure, target_structure, radius=2.0)
    final_score = descriptor_score
    if ca_clashes > ca_clash_threshold or heavy_clashes > heavy_atom_clash_threshold:
        final_score = 0.0
    return ca_clashes, heavy_clashes, final_score


def compute_score_and_clashes(
    P1_pdb,
    P1_pcd,
    P1_descs,
    P1_site,
    P1_indices,
    P2_pdb,
    P2_pcd,
    P2_descs,
    P2_site,
    P2_indices,
    compute_clashes: bool,
    **kwargs,
):
    """Descriptor score and optional clash counts between two aligned patches."""
    P1_patch_coords = np.array(P1_pcd.points)[P1_indices[P1_site]]
    normalized_desc_dist_score = compute_descriptor_score(
        P1_patch_coords,
        P1_descs,
        P1_indices,
        P1_site,
        P2_pcd,
        P2_descs,
        P2_indices,
        P2_site,
    )
    output = [[0, 0], normalized_desc_dist_score]

    if compute_clashes:
        assert "target_structure" in kwargs, "Please provide target_structure for clash counting."
        assert "target_chain" in kwargs, "Please provide target_chain for clash counting."
        target_structure = get_filtered_target_structure(
            kwargs["target_structure"],
            kwargs["target_chain"],
            kwargs.get("_target_cache", {}),
        )
        ca_clashes, heavy_clashes, final_score = compute_hit_clash_score(
            P1_pdb,
            target_structure,
            normalized_desc_dist_score,
            ca_clash_threshold=kwargs.get("ca_clash_threshold", 1.0),
            heavy_atom_clash_threshold=kwargs.get("heavy_atom_clash_threshold", 5.0),
        )
        output[0][0], output[0][1] = ca_clashes, heavy_clashes
        output[1] = final_score
        return output, target_structure
    return output, None

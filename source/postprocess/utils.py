import math
import os
from pathlib import Path
from io import StringIO
from typing import Union
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.SVDSuperimposer import SVDSuperimposer
from sklearn.cluster import DBSCAN
from tqdm import tqdm

from src.metrics.secondary_structure import find_sse


DATABASE_ROOT = os.environ.get(
    "DATABASE_ROOT",
    os.environ.get(
        "DATA_ROOT",
        "/work/lpdi/users/diazrovi/domaindome/20260221-AFDBv6_domaindome_DPAM_masif/dpam_domaindome_masif_db",
    ),
)

# 20 standard amino acid three-letter codes
STANDARD_AA = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
})


def filter_structure_to_standard_aa(struct, standard_aa=None):
    """
    Return a new Structure with only standard amino acid residues.
    Preserves chain IDs and residue numbering.
    """
    if standard_aa is None:
        standard_aa = STANDARD_AA

    model = struct[0] if isinstance(struct, Structure) else struct
    builder = StructureBuilder()
    builder.init_structure('filtered')
    builder.init_model(0)

    for chain in model.get_chains():
        builder.init_chain(chain.id)

    out_struct = builder.get_structure()
    out_model = out_struct[0]

    for chain in model.get_chains():
        chain_id = chain.id
        for res in chain.get_residues():
            if res.get_resname() in standard_aa:
                out_model[chain_id].add(res.copy())

    return out_struct


def maybe_load_structure(path_or_structure: Union[str, Path, Structure]):
    return path_or_structure if isinstance(path_or_structure, Structure) else PDBParser(QUIET=True).get_structure('', str(path_or_structure))


def get_descriptor(row, descriptor_root=Path(DATABASE_ROOT) / 'descriptors' / 'sc05' / 'all_feat', flipped=False):
    if math.isnan(row.matched_vix):
        return None
    desc_file = Path(descriptor_root, row.matched_protein, 'p1_desc_flipped.npy' if flipped else 'p1_desc_straight.npy')
    return np.load(desc_file)[int(row.matched_vix)]


def get_transformed_struct(row, database_root=DATABASE_ROOT):
    original_pdb = Path(database_root, 'data_preparation', '01-benchmark_pdbs', row.matched_protein + '.pdb')
    struct = PDBParser().get_structure(row.matched_protein, original_pdb)

    transform = np.array(list(map(float, row.flattened_transform.split(',')))).reshape(4, 4)
    struct.transform(rot=transform[:3, :3].T, tran=transform[:3, 3])
    return struct


def get_pdb(row, database_root=DATABASE_ROOT):

    struct = get_transformed_struct(row, database_root)

    out = StringIO()
    pdb_io = PDBIO()
    pdb_io.set_structure(struct)
    pdb_io.save(out)

    return out.getvalue()


class GLoopMatcher:
    # BB_ATOMS = ['N', 'CA', 'C', 'O']
    BB_ATOMS = ['N', 'CA', 'C']

    def __init__(self, template_chain, template_gly, rel_resi_range=(-5, 2), max_gly_dist=4.0):

        glycine = template_chain[template_gly]
        assert glycine.get_resname() == 'GLY'
        self.gly_pos = glycine['CA'].get_coord().reshape(1, 3)
        self.max_gly_dist = max_gly_dist
        self.rel_start = rel_resi_range[0]
        self.rel_stop = rel_resi_range[1]

        start, stop = template_gly + self.rel_start, template_gly + self.rel_stop
        assert start < template_gly < stop
        self.template_backbone = np.stack(
            [template_chain[resi][a].get_coord() for resi in range(start, stop + 1) for a in self.BB_ATOMS]
        )

    def __call__(self, struct):
        residues = {res.id[1]: res for res in struct.get_residues()}

        # Look for glycine near reference position
        gly_residues = [r for r in struct.get_residues() if r.get_resname() == 'GLY']
        if len(gly_residues) < 1:
            return None  # no glycine
        ca_coords = np.stack([res['CA'].get_coord() for res in gly_residues])
        dists = np.linalg.norm(ca_coords - self.gly_pos, axis=-1)
        candidates = np.where(dists <= self.max_gly_dist)[0]
        if len(candidates) < 1:
            return None  # no glycine close enough

        rmsd_vals = []
        for cidx in candidates:
            gly_resi = gly_residues[cidx].id[1]
            try:
                motif = [residues[gly_resi + offset] for offset in range(self.rel_start, self.rel_stop + 1)]
            except KeyError as e:
                # print('one or more residues not found', e)
                continue  # one or more residues not found

            if not all(motif[i].id[1] == motif[i - 1].id[1] + 1 for i in range(1, len(motif))):
                # print('sequence gap')
                continue  # sequence gap

            # Align local peptide stretch to template G-loop
            try:
                backbone = np.stack([res[a].get_coord() for res in motif for a in self.BB_ATOMS])
            except KeyError as e:
                # print("Atom not found", e)
                continue
            assert len(backbone) == len(self.template_backbone)
            sup = SVDSuperimposer()
            sup.set(self.template_backbone, backbone)
            sup.run()

            # Compute RMSD
            rmsd_vals.append(sup.get_rms())

        return min(rmsd_vals) if len(rmsd_vals) > 0 else None


def structural_clusters(df, rmsd_thresh=5.0, database_root=DATABASE_ROOT):
    """
    Clusters domains according to their binding mode (RMSD-based).
    Adds the following new columns to the input dataframe:
    - cluster_id: unique cluster id
    - cluster_size: number members in the same cluster
    - cluster_mean_rmsd: mean RMSD to all other cluster members
    """

    def get_coords(struct, atoms_to_keep=['N', 'CA', 'C']):
        return np.stack([a.get_coord() for a in struct.get_atoms() if a.name in atoms_to_keep])


    def rmsd(coords1, coords2):
        assert coords1.shape == coords2.shape
        diff = coords1 - coords2
        return np.sqrt(np.mean(np.sum(diff * diff, axis=-1)))


    current_cluster_label = 0
    # for domain_name in ['Q5IJ48-F1-dom-01_A']:
    for domain_name in tqdm(df.matched_protein.unique()):

        domain_table = df[df.matched_protein == domain_name]

        if len(domain_table) == 1:
            df.loc[domain_table.index, 'cluster_id'] = current_cluster_label
            df.loc[domain_table.index, 'cluster_size'] = 1
            df.loc[domain_table.index, 'cluster_mean_rmsd'] = 0.0
            current_cluster_label += 1
            continue

        positions = [get_coords(get_transformed_struct(row, database_root=database_root)) for i, row in domain_table.iterrows()]

        rmsd_vals = np.zeros((len(domain_table), len(domain_table)))
        for i in range(len(domain_table)):
            for j in range(i + 1, len(domain_table)):
                rmsd_ij = rmsd(positions[i], positions[j])
                rmsd_vals[i, j] = rmsd_ij
                rmsd_vals[j, i] = rmsd_ij

        # DBSCAN clustering
        labels = DBSCAN(eps=rmsd_thresh, min_samples=2, metric='precomputed').fit_predict(rmsd_vals)

        # Create new clusters for outliers
        labels[labels == -1] = np.arange((labels == -1).sum()) + labels.max() + 1

        for lb in set(labels):
            cluster_table = domain_table.iloc[labels == lb]
            df.loc[cluster_table.index, 'cluster_id'] = current_cluster_label
            df.loc[cluster_table.index, 'cluster_size'] = len(cluster_table)
            df.loc[cluster_table.index, 'cluster_mean_rmsd'] = rmsd_vals[labels == lb][:, labels == lb].mean(axis=1)
            current_cluster_label += 1

    return df


# def get_location_masks(row):

#     extracellular = np.array([x == '1' for x in row.extracellular_mask]) if isinstance(row.extracellular_mask, str) else None
#     transmembrane = np.array([x == '0' for x in row.non_transmembrane_mask]) if isinstance(row.non_transmembrane_mask, str) else None

#     if extracellular is not None and transmembrane is not None:
#         intracellular = ~transmembrane & ~extracellular
#     else:
#         intracellular = None

#     return intracellular, extracellular, transmembrane
def get_location_masks(row):

    if not isinstance(row.deeptmhmm_annotation, str):
        return None, None, None, None

    # signal peptide (S) or inside cell/cytosol (I)
    intracellular = np.array([label in {'S', 'I'} for label in row.deeptmhmm_annotation])

    # outside cell/lumen of ER/Golgi/lysosomes (O)
    extracellular = np.array([label in {'O'} for label in row.deeptmhmm_annotation])

    # alpha membrane (M) or beta membrane (B)
    transmembrane = np.array([label in {'M', 'B'} for label in row.deeptmhmm_annotation])

    # periplasm (P) or (?)
    unassigned = np.array([label in {'P', '?'} for label in row.deeptmhmm_annotation])

    return intracellular, extracellular, transmembrane, unassigned


def iface_sse_metrics(row, iface_resi):
    """
    Calculate secondary structure element metrics for interface residues.
    [Note: exclude 3 residues from each terminus of mdtraj_dssp_annotation]
    
    Args:
        row: DataFrame row containing mdtraj_dssp_annotation and resi columns
        iface_resi: Set of interface residue numbers
    
    Returns:
        Tuple of (total_n_sse, n_helix, n_strand, helix_frac, strand_frac)
    """
    if not isinstance(row.mdtraj_dssp_annotation, str):
        return None, None, None, None, None
    
    # Remove the first and last 3 letters of the mdtraj_dssp_annotation, before counting the number of secondary structure elements
    mdtraj_dssp_annotation = row.mdtraj_dssp_annotation[3:-3]

    residue_ids = list(map(int, row.resi.split(',')))
    mdtraj_dssp_labels = list(row.mdtraj_dssp_annotation)
    assert len(residue_ids) == len(mdtraj_dssp_labels)

    iface_mdtraj_dssp_labels = [x for x, resi in zip(mdtraj_dssp_labels, residue_ids) if resi in iface_resi]

    segments = find_sse(residue_ids, mdtraj_dssp_labels)
    iface_segments = [seg for seg in segments if any(seg['start'] <= resi <= seg['end'] for resi in iface_resi)]

    total_n_sse = len(iface_segments)
    n_helix = sum(s['label'] == 'H' for s in iface_segments)
    n_strand = sum(s['label'] == 'E' for s in iface_segments)

    helix_frac = np.mean([x == 'H' for x in iface_mdtraj_dssp_labels])
    strand_frac = np.mean([x == 'E' for x in iface_mdtraj_dssp_labels])

    return total_n_sse, n_helix, n_strand, helix_frac, strand_frac


def iface_metrics(row, iface_col='matched_iface_resi'):
    """
    Aggregate metrics for interface residues.
    
    Args:
        row: DataFrame row containing matched_iface_resi (or custom interface column), 
             resi, plddt, deeptmhmm_annotation, and mdtraj_dssp_annotation columns
        iface_col: Name of the column containing interface residues (default: 'matched_iface_resi')
    
    Returns:
        Tuple of interface metrics: (iface_intracellular_frac, iface_extracellular_frac,
        iface_transmembrane_frac, iface_plddt, iface_n_sse, iface_n_helix, iface_n_strand,
        iface_helix_frac, iface_strand_frac)
    """
    resi = np.array(list(map(int, row.resi.split(','))))

    if type(getattr(row, iface_col)) != str:
        # if no interface residues were found
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    iface_resi = set(map(int, getattr(row, iface_col).split(',')))
    iface_mask = np.array([r in iface_resi for r in resi])

    # Average pLDDT at the interface
    plddt = np.array(list(map(float, row.plddt.split(','))))
    assert len(plddt) == len(resi)
    iface_plddt = np.mean(plddt[iface_mask])

    # Fraction of residues labelled as extracellular at the interface
    intracellular_mask, extracellular_mask, transmembrane_mask, _ = get_location_masks(row)

    if intracellular_mask is not None:
        assert len(resi) == len(intracellular_mask)
        intracellular_resi = set(resi[intracellular_mask])
        iface_intracellular_frac = len(iface_resi & intracellular_resi) / len(iface_resi)
    else:
        iface_intracellular_frac = None

    if extracellular_mask is not None:
        assert len(resi) == len(extracellular_mask)
        extracellular_resi = set(resi[extracellular_mask])
        iface_extracellular_frac = len(iface_resi & extracellular_resi) / len(iface_resi)
    else:
        iface_extracellular_frac = None

    # Fraction of residues labelled as non-transmembrane at the interface
    if transmembrane_mask is not None:
        assert len(resi) == len(transmembrane_mask)
        transmembrane_resi = set(resi[transmembrane_mask])
        iface_transmembrane_frac = len(iface_resi & transmembrane_resi) / len(iface_resi)
    else:
        iface_transmembrane_frac = None

    # Secondary structure elements at the interface
    iface_n_sse, iface_n_helix, iface_n_strand, iface_helix_frac, iface_strand_frac = iface_sse_metrics(row, iface_resi)

    return iface_intracellular_frac, iface_extracellular_frac, iface_transmembrane_frac, iface_plddt, iface_n_sse, iface_n_helix, iface_n_strand, iface_helix_frac, iface_strand_frac


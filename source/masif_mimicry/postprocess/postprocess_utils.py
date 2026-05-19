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

from masif_mimicry.config.paths import resolve_database_paths
from masif_mimicry.postprocess.metrics.secondary_structure import find_sse
from masif_mimicry.postprocess.structures import (
    STANDARD_AA,
    filter_structure_to_standard_aa,
    maybe_load_structure,
)
from masif_mimicry.utils.transforms import get_transformed_struct_from_row


def _descriptor_root(database_dir):
    db_root, _ = resolve_database_paths(database_dir)
    return Path(db_root) / "descriptors" / "sc05" / "all_feat"


def get_descriptor(row, database_dir, flipped=False):
    if math.isnan(row.matched_vix):
        return None
    desc_file = Path(
        _descriptor_root(database_dir),
        row.matched_protein,
        "p1_desc_flipped.npy" if flipped else "p1_desc_straight.npy",
    )
    return np.load(desc_file)[int(row.matched_vix)]


def get_transformed_struct(row, database_dir):
    """Apply flattened_transform to a domainome PDB (legacy matched_protein rows)."""
    if hasattr(row, "P1_id") and hasattr(row, "flattened_transform"):
        return get_transformed_struct_from_row(row, database_dir)

    _, db_prep = resolve_database_paths(database_dir)
    original_pdb = Path(db_prep, "01-benchmark_pdbs", f"{row.matched_protein}.pdb")
    struct = PDBParser().get_structure(row.matched_protein, original_pdb)
    transform = np.array(list(map(float, row.flattened_transform.split(",")))).reshape(4, 4)
    struct.transform(rot=transform[:3, :3].T, tran=transform[:3, 3])
    return struct


def get_pdb(row, database_dir):
    struct = get_transformed_struct(row, database_dir)

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


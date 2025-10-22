#!/usr/bin/env python3
"""
proc_trunc_masif_hits.py

This script processes a CSV table (created by postprocess_masif_hits.py) to compute structural and interface
metrics for truncated binder PDB files. For each row, the target PDB (from target_path) and the truncated binder
PDB (from trunc_binder_path) are loaded and analyzed. The target.vert file is generated using the matched_protein_path.
Metrics computed for the truncated binder are renamed according to the mapping below:

    Original                      Truncated
    ---------------------------------------------------
    ligand_binder_buried_SASA  ->  ligand_trunc_buried_SASA
    binder_iface_n_resi       ->  trunc_iface_n_resi
    binder_iface_plddt        ->  trunc_iface_plddt
    binder_iface_residues     ->  trunc_iface_residues
    binder_iface_intracellular->  trunc_iface_intracellular
    binder_resi_start         ->  trunc_resi_start
    binder_resi_end           ->  trunc_resi_end
    binder_length             ->  trunc_length
    binder_total_n_ss         ->  trunc_total_n_ss
    binder_n_structured_ss    ->  trunc_n_structured_ss
    binder_structured_percent ->  trunc_structured_percent
    binder_iface_n_ss         ->  trunc_iface_n_ss
    binder_iface_structured   ->  trunc_iface_structured
    binder_iface_ss           ->  trunc_iface_ss
    abs_geodesic_length       ->  abs_geodesic_length_trunc
    norm_geodesic_length      ->  norm_geodesic_length_trunc

In addition, new columns are added as differences between the truncated and original binder interface metrics:
    - d_iface_n_resi
    - d_iface_n_ss
    - d_iface_structured

A new metric, EvoEF2_score, is also extracted from the trunc_binder_path filename. For example, given:
    /scratch/ymeng/MaSIF_search_3.6/8VLB_A_A_3JF_301/Truncate_86/out/Q8NA31-F1-dom-01_A_1465_101_186_-451.pdb
the value "-451" (converted to a numeric value) is assigned to EvoEF2_score.

Usage:
    python proc_trunc_masif_hits.py --input input.csv --out_csv_file output.csv [--ligand CHAIN_RESNAME] [--sulfur ATOM_NAME] [--debug]
"""

from pathlib import Path
from argparse import ArgumentParser
import subprocess
import math
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.DSSP import DSSP
from Bio.PDB.Polypeptide import protein_letters_3to1  # for sequence extraction
from Bio import pairwise2  # for alignment
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB import StructureBuilder
import yaml
import sys
# Add the script directory to PATH to import geodesic_length
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(script_dir)
from geodesic_length import compute_geodesic_length

# ------------------ Load configuration ------------------
# Load configuration from config/config.yaml
CONFIG_PATH = os.path.join(os.path.dirname(__file__), "config.yaml")
with open(CONFIG_PATH, 'r') as f:
    config = yaml.safe_load(f)

paths_config = config['paths']
deeptmhmm_dir = Path(paths_config['deeptmhmm_dir'])
stride_exec = paths_config['stride_exec']


# ------------------ Helper Functions ------------------

def coord_and_radii(rdmol, ignore={'H'}):
    _periodic_table = Chem.GetPeriodicTable()
    coord = rdmol.GetConformer().GetPositions()
    radii = np.array([_periodic_table.GetRvdw(a.GetSymbol()) for a in rdmol.GetAtoms()])
    mask = np.array([a.GetSymbol() not in ignore for a in rdmol.GetAtoms()])
    coord = coord[mask]
    radii = radii[mask]
    assert coord.shape[0] == radii.shape[0]
    return coord, radii

def count_clashes(pdb1, pdb2, strictness=1.0):
    """Compute number of clashing heavy atom pairs."""
    protein1 = Chem.MolFromPDBFile(str(pdb1), sanitize=False)
    protein2 = Chem.MolFromPDBFile(str(pdb2), sanitize=False)
    coord1, radii1 = coord_and_radii(protein1)
    coord2, radii2 = coord_and_radii(protein2)
    dist = np.sqrt(np.sum((coord1[:, None, :] - coord2[None, :, :]) ** 2, axis=-1))
    clashes = dist < strictness * (radii1[:, None] + radii2[None, :])
    num_clashes = np.sum(clashes)
    return num_clashes

def extract_ligand_residue(pdb_struct, ligand_chain, ligand_name):
    ligand = [res for res in pdb_struct[ligand_chain].get_residues() if res.get_resname() == ligand_name]
    assert len(ligand) == 1
    return ligand[0]

def copy_residues(input_residues, output_struct, output_chain):
    if not isinstance(input_residues, list):
        input_residues = list(input_residues.get_residues())
    for i, res in enumerate(input_residues):
        _res = res.copy()
        _res.id = ("", i + 1, "")
        output_struct[output_chain].add(_res)

def merge_structures(target_struct, binder_struct, ligand=None):
    builder = StructureBuilder.StructureBuilder()
    builder.init_structure('complex')
    builder.init_model(0)
    builder.init_chain("T")  # for target
    builder.init_chain("B")  # for binder
    if ligand is not None:
        builder.init_chain("L")  # for ligand
    out_struct = builder.get_structure()[0]
    copy_residues(target_struct, out_struct, output_chain="T")
    copy_residues(binder_struct, out_struct, output_chain="B")
    if ligand is not None:
        ligand_residue = extract_ligand_residue(target_struct, ligand["chain"], ligand["name"])
        copy_residues([ligand_residue], out_struct, output_chain="L")
    return out_struct

def merge_binder_ligand(binder_struct, ligand_residue):
    """Builds a structure containing only the binder and the ligand."""
    builder = StructureBuilder.StructureBuilder()
    builder.init_structure('binder_ligand_complex')
    builder.init_model(0)
    builder.init_chain("B")
    builder.init_chain("L")
    new_struct = builder.get_structure()[0]
    copy_residues(binder_struct, new_struct, output_chain="B")
    copy_residues([ligand_residue], new_struct, output_chain="L")
    return new_struct

def compute_sasa_values(target_pdb, binder_pdb, ligand_def: dict = None):
    out = {}
    target = PDBParser(QUIET=True).get_structure('', target_pdb)[0]
    binder = PDBParser(QUIET=True).get_structure('', binder_pdb)[0]
    whole_complex = merge_structures(target, binder, ligand_def)

    # Compute SASA for target
    ShrakeRupley().compute(target, level="R")
    out['target_sasa'] = sum(r.sasa for r in target.get_residues())

    # Binder SASA
    ShrakeRupley().compute(binder, level="R")
    out['binder_sasa'] = sum(r.sasa for r in binder.get_residues())

    # Complex SASA
    ShrakeRupley().compute(whole_complex, level="R")
    out['complex_sasa'] = sum(r.sasa for r in whole_complex.get_residues())
    target_residues_in_complex = list(whole_complex["T"].get_residues())
    out['target_complex_sasa'] = sum(r.sasa for r in target_residues_in_complex)
    out['target_buried_sasa'] = out['target_sasa'] - out['target_complex_sasa']
    binder_residues_in_complex = list(whole_complex["B"].get_residues())
    out['binder_complex_sasa'] = sum(r.sasa for r in binder_residues_in_complex)
    out['binder_buried_sasa'] = out['binder_sasa'] - out['binder_complex_sasa']

    if ligand_def is not None:
        ligand = extract_ligand_residue(target, ligand_def["chain"], ligand_def["name"])
        ligand_in_complex = extract_ligand_residue(whole_complex, "L", ligand_def["name"])

        def compute_isolated_residue_sasa(res):
            builder = StructureBuilder.StructureBuilder()
            builder.init_structure('isolated')
            builder.init_model(0)
            builder.init_chain("I")
            new_struct = builder.get_structure()[0]
            res_copy = res.copy()
            res_copy.id = ("", 1, "")
            new_struct["I"].add(res_copy)
            ShrakeRupley().compute(new_struct, level="R")
            return res_copy.sasa

        ligand_isolated_sasa = compute_isolated_residue_sasa(ligand)
        ligand_complex_sasa = ligand_in_complex.sasa

        out['ligand_sasa'] = ligand.sasa
        out['ligand_complex_sasa'] = ligand_complex_sasa
        out['ligand_buried_sasa'] = out['ligand_sasa'] - ligand_complex_sasa

        binder_ligand_complex = merge_binder_ligand(binder, ligand)
        ShrakeRupley().compute(binder_ligand_complex, level="R")
        ligand_in_binder_complex = extract_ligand_residue(binder_ligand_complex, "L", ligand_def["name"])
        binder_ligand_complex_sasa = ligand_in_binder_complex.sasa

        out['ligand_binder_buried_SASA'] = ligand_isolated_sasa - binder_ligand_complex_sasa
        out['ligand_iface_sasa_contribution'] = out['ligand_buried_sasa'] / out['target_buried_sasa']

    return out

def get_contact_atoms(target_struct, ligand_def=None):
    """Extract heavy atoms from target (and ligand, if defined) for interface contact."""
    contact_atoms = []
    if ligand_def is not None:
        ligand_res = extract_ligand_residue(target_struct, ligand_def["chain"], ligand_def["name"])
        contact_atoms.extend([atom for atom in ligand_res if atom.element != 'H'])
        for chain in target_struct:
            for res in chain:
                if res != ligand_res:
                    for atom in res:
                        if atom.element != 'H':
                            contact_atoms.append(atom)
    else:
        contact_atoms = [atom for atom in target_struct.get_atoms() if atom.element != 'H']
    return contact_atoms


def get_interface_residues(binder_struct, contact_atoms, cutoff=4.5):
    """Return binder residues that have any heavy atom within the cutoff distance."""
    ns = NeighborSearch(contact_atoms)
    iface_residues = []
    for res in binder_struct.get_residues():
        for atom in res:
            if atom.element == 'H':
                continue
            if ns.search(atom.get_coord(), cutoff):
                iface_residues.append(res)
                break
    return iface_residues

def compute_iface_metrics(binder_struct, iface_residues):
    """Compute average CA B-factor for interface residues."""
    total_bfactor = 0.0
    count_atoms = 0
    for res in iface_residues:
        for atom in res:
            if atom.element == 'H' or atom.name != "CA":
                continue
            total_bfactor += atom.get_bfactor()
            count_atoms += 1
    avg_bfactor = total_bfactor / count_atoms if count_atoms > 0 else float('nan')
    return avg_bfactor

def filter_contact_atoms(contact_atoms, filter_by_atoms, cutoff=6.0):
    """
    Filter the list of contact atoms to only those within the specified cutoff distance
    from any atom in filter_by_atoms.

    Args:
        contact_atoms (list): List of Bio.PDB Atom objects to be filtered.
        filter_by_atoms (list): List of Bio.PDB Atom objects used for distance filtering.
        cutoff (float): Distance cutoff in Angstroms. Default is 6.0.

    Returns:
        list: Filtered list of atoms.
    """
    import numpy as np

    # Precompute coordinates of the atoms used for filtering.
    filter_coords = np.array([atom.get_coord() for atom in filter_by_atoms])
    filtered_atoms = []
    for atom in contact_atoms:
        atom_coord = atom.get_coord()  # shape: (3,)
        distances = np.linalg.norm(filter_coords - atom_coord, axis=1)
        if np.any(distances < cutoff):
            filtered_atoms.append(atom)
    return filtered_atoms

def compute_binder_interface_metrics(binder_struct, target_struct, binder_pdb_path, ligand_def, debug=False):
    """
    Compute binder interface metrics including:
      - Number of interface residues and average CA B-factor.
      - Interface residue numbering.
      - Intracellular interface prediction using DeepTMHMM.
      - Residue range and STRIDE-based secondary structure metrics.
      - Secondary structure metrics for the binder interface.
    
    This modified version filters the contact atoms to those within 6 Å of any ligand atom.
    
    Returns a dictionary with the computed metrics.
    """
    metrics = {}
    try:
        # Primary computation: get contact atoms from target structure.
        contact_atoms = get_contact_atoms(target_struct, ligand_def)
        
        # If a ligand is defined, filter contact atoms to those near the ligand atoms.
        if ligand_def is not None:
            # Extract ligand residues and then get its non-hydrogen atoms.
            ligand_res = extract_ligand_residue(target_struct, ligand_def["chain"], ligand_def["name"])
            ligand_atoms = [atom for atom in ligand_res if atom.element != 'H']
            # Use the modified filtering function with a cutoff of 6 Å.
            contact_atoms = filter_contact_atoms(contact_atoms, ligand_atoms, cutoff=6.0)
        
        # Compute interface residues for the binder structure.
        iface_residues = get_interface_residues(binder_struct, contact_atoms, cutoff=4.5)
        metrics['binder_iface_n_resi'] = len(iface_residues)
        avg_bfactor = compute_iface_metrics(binder_struct, iface_residues)
        metrics['binder_iface_plddt'] = avg_bfactor
        metrics['binder_iface_residues'] = '_'.join(str(res.id[1]) for res in iface_residues)

        # Compute binder_iface_intracellular using DeepTMHMM prediction.
        try:
            binder_pdb_file = binder_pdb_path
            uniprot_id = binder_pdb_file.name.split("-")[0]
            deeptmhmm_path = deeptmhmm_dir / uniprot_id / "predicted_topologies.3line"
            if not deeptmhmm_path.exists():
                if debug:
                    print(f"No .3line file found for binder {binder_pdb_file.name}")
                binder_iface_intracellular = 1
            else:
                with open(deeptmhmm_path, "r") as f:
                    lines = f.readlines()
                if len(lines) < 3:
                    raise ValueError("DeepTMHMM file format incorrect.")
                dt_seq = lines[1].strip()
                dt_pred = lines[2].strip()

                binder_seq = ""
                residue_mapping = {}
                seq_index = 0
                for chain in binder_struct:
                    for res in chain:
                        if res.id[0] != " ":
                            continue
                        try:
                            aa = protein_letters_3to1(res.get_resname())
                        except Exception:
                            aa = "-"
                            continue
                        binder_seq += aa
                        residue_mapping[(chain.id, res.id)] = seq_index
                        seq_index += 1

                alignments = pairwise2.align.globalms(binder_seq, dt_seq, 2, -1, -0.5, -0.1)
                best_alignment = alignments[0]
                aligned_binder = best_alignment[0]
                aligned_dt = best_alignment[1]
                dt_pred_aligned = ""
                dt_index = 0
                for a in aligned_dt:
                    if a == "-":
                        dt_pred_aligned += "-"
                    else:
                        dt_pred_aligned += dt_pred[dt_index]
                        dt_index += 1

                n_non_gap = sum(1 for a in aligned_binder if a != "-")
                n_matches = sum(
                    1 for i in range(len(aligned_binder))
                    if aligned_binder[i] != "-" and aligned_binder[i] == aligned_dt[i]
                )
                identity = (n_matches / n_non_gap) * 100 if n_non_gap > 0 else 0

                if debug:
                    print(f"Alignment for {binder_pdb_file.name}:")
                    print(aligned_binder)
                    print(aligned_dt)
                    print(dt_pred_aligned)
                    print(f"Identity: {identity:.2f}%")
                
                if identity < 95:
                    binder_iface_intracellular = 1
                else:
                    unaligned_to_aligned = {}
                    unaligned_idx = 0
                    for aligned_idx, char in enumerate(aligned_binder):
                        if char != "-":
                            unaligned_to_aligned[unaligned_idx] = aligned_idx
                            unaligned_idx += 1

                    intracellular_count = 0
                    total_count = 0
                    for res in iface_residues:
                        key = (res.get_parent().id, res.id)
                        if key not in residue_mapping:
                            continue
                        binder_unaligned_index = residue_mapping[key]
                        if binder_unaligned_index not in unaligned_to_aligned:
                            continue
                        aligned_index = unaligned_to_aligned[binder_unaligned_index]
                        if aligned_index < len(aligned_dt) and aligned_binder[aligned_index] == aligned_dt[aligned_index]:
                            pred = dt_pred_aligned[aligned_index]
                        else:
                            pred = "I"
                        total_count += 1
                        if pred not in ["O", "M"]:
                            intracellular_count += 1
                    binder_iface_intracellular = intracellular_count / total_count if total_count > 0 else 1
            metrics['binder_iface_intracellular'] = binder_iface_intracellular
        except Exception as e:
            binder_name = binder_pdb_file.name if 'binder_pdb_file' in locals() else ""
            if debug:
                print(f"Error computing binder_iface_intracellular (primary) for {binder_name}: {e}")
            metrics['binder_iface_intracellular'] = 1

        # Compute binder residue range metrics.
        binder_res_nums = []
        for chain in binder_struct:
            for res in chain:
                if res.id[0] != " ":
                    continue
                binder_res_nums.append(res.id[1])
        if binder_res_nums:
            metrics['binder_resi_start'] = min(binder_res_nums)
            metrics['binder_resi_end'] = max(binder_res_nums)
            metrics['binder_length'] = metrics['binder_resi_end'] - metrics['binder_resi_start'] + 1
        else:
            metrics['binder_resi_start'] = None
            metrics['binder_resi_end'] = None
            metrics['binder_length'] = None

        # STRIDE-based secondary structure metrics.
        if metrics['binder_resi_start'] is not None:
            try:
                binder_pdb_file = binder_pdb_path
                result = subprocess.run([stride_exec, str(binder_pdb_file)], capture_output=True, text=True, check=True)
                stride_output = result.stdout

                # Build a mapping from residue number to secondary structure assignment.
                stride_mapping = {}
                for line in stride_output.splitlines():
                    if line.startswith("ASG"):
                        parts = line.split()
                        try:
                            res_num = int(parts[3])
                        except ValueError:
                            continue
                        ss_code = parts[5]
                        ss_code = ss_code if ss_code in ['H', 'E'] else 'L'
                        stride_mapping[res_num] = ss_code

                all_residues = []
                for chain in binder_struct:
                    for res in chain:
                        if res.id[0] != " ":
                            continue
                        resi = res.id[1]
                        resn = res.get_resname()
                        ss_type = stride_mapping.get(resi, "L")
                        all_residues.append((resi, resn, ss_type))
                all_residues.sort(key=lambda x: x[0])

                seg_data = []
                if all_residues:
                    current_seg = 1
                    seg_data.append((all_residues[0][0], all_residues[0][1], all_residues[0][2], current_seg))
                    for i in range(1, len(all_residues)):
                        prev = all_residues[i-1]
                        curr = all_residues[i]
                        if curr[0] == prev[0] + 1 and curr[2] == prev[2]:
                            seg_data.append((curr[0], curr[1], curr[2], current_seg))
                        else:
                            current_seg += 1
                            seg_data.append((curr[0], curr[1], curr[2], current_seg))
                else:
                    current_seg = 0

                unique_segments = {}
                for seg in seg_data:
                    seg_index = seg[3]
                    if seg_index not in unique_segments:
                        unique_segments[seg_index] = seg[2]
                binder_total_n_ss = len(unique_segments)
                binder_n_structured_ss = sum(1 for s in unique_segments.values() if s in ['H', 'E'])
                total_residues = len(all_residues)
                structured_residues_list = [res for res in all_residues if res[2] in ['H', 'E']]
                structured_residues = len(structured_residues_list)
                binder_structured_percent = (structured_residues / total_residues * 100) if total_residues > 0 else 0
                metrics['binder_total_n_ss'] = binder_total_n_ss
                metrics['binder_n_structured_ss'] = binder_n_structured_ss
                metrics['binder_structured_percent'] = binder_structured_percent

                interface_res_nums = set(res.id[1] for res in iface_residues)
                iface_seg_data = [row for row in seg_data if row[0] in interface_res_nums]
                iface_seg_indices = {row[3] for row in iface_seg_data}
                n_ss_iface = len(iface_seg_indices)
                total_iface = len(iface_seg_data)
                structured_iface = sum(1 for row in iface_seg_data if row[2] in ['H', 'E'])
                structured_fraction = structured_iface / total_iface if total_iface > 0 else 0

                if debug:
                    print("resi\tresn\tss_type\tss_index\tiface")
                    for row in seg_data:
                        iface_marker = "*" if row[0] in interface_res_nums else ""
                        print(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{iface_marker}")
                    print(f"n_ss: {n_ss_iface}")
                    print(f"structured: {structured_fraction:.2f}")

                metrics['binder_iface_n_ss'] = n_ss_iface
                metrics['binder_iface_structured'] = structured_fraction

                ss_map = {resi: ss for resi, resn, ss in all_residues}
                binder_iface_ss = "_".join(ss_map.get(res.id[1], 'L') for res in iface_residues if res.id[0] == " ")
                metrics['binder_iface_ss'] = binder_iface_ss
            except Exception as e:
                if debug:
                    print(f"Error computing STRIDE-based interface metrics: {e}")
                metrics['binder_iface_n_ss'] = None
                metrics['binder_iface_structured'] = None
                metrics['binder_iface_ss'] = None
        else:
            metrics['binder_iface_n_ss'] = None
            metrics['binder_iface_structured'] = None
            metrics['binder_iface_ss'] = None

    except Exception as primary_error:
        print(f"Warning: Primary binder interface metrics computation failed for {binder_pdb_path}: {primary_error}")
        metrics['binder_iface_n_resi'] = None
        metrics['binder_iface_plddt'] = None     
        metrics['binder_iface_intracellular'] = 1
        metrics['binder_resi_start'] = None
        metrics['binder_resi_end'] = None
        metrics['binder_length'] = None
        metrics['binder_iface_n_ss'] = None
        metrics['binder_iface_structured'] = None
        metrics['binder_iface_ss'] = None

    return metrics




# ------------------ Processing Function for Truncated Binder Metrics ------------------

def process_truncated_results(input_csv, ligand_def, sulfur_name, debug):
    df = pd.read_csv(input_csv)
    results = []
    parser = PDBParser(QUIET=True)
    
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing CSV rows"):
        # Get paths from CSV (ensuring they are absolute)
        try:
            target_path = Path(row['target_path']).resolve()
            trunc_binder_path = Path(row['trunc_binder_path']).resolve()
            matched_protein_path = Path(row['matched_protein_path']).resolve()
        except Exception as e:
            print(f"Error reading file paths in row {idx}: {e}")
            continue

        # Parse structures
        try:
            target_struct = parser.get_structure("", str(target_path))[0]
        except Exception as e:
            print(f"Error parsing target structure {target_path}: {e}")
            continue
        try:
            binder_struct = parser.get_structure("", str(trunc_binder_path))[0]
        except Exception as e:
            print(f"Error parsing truncated binder structure {trunc_binder_path}: {e}")
            continue

        # Recompute clash metrics using truncated binder
        clashes1 = count_clashes(str(target_path), str(trunc_binder_path), strictness=1.0)
        clashes075 = count_clashes(str(target_path), str(trunc_binder_path), strictness=0.75)

        # Recompute SASA metrics
        sasa_metrics = compute_sasa_values(str(target_path), str(trunc_binder_path), ligand_def=ligand_def)
        # Add new metric: trunc_sasa (SASA for the truncated binder structure)
        sasa_metrics["trunc_sasa"] = sasa_metrics.get("binder_sasa", None)
        # Rename ligand binder metric to truncated version if present.
        if "ligand_binder_buried_SASA" in sasa_metrics:
            sasa_metrics["ligand_trunc_buried_SASA"] = sasa_metrics.pop("ligand_binder_buried_SASA")

        # Recompute binder interface metrics on truncated binder.
        binder_metrics_trunc = compute_binder_interface_metrics(binder_struct, target_struct, trunc_binder_path, ligand_def, debug=debug)

        # Rename binder interface metrics to new truncated column names.
        mapping = {
            'binder_iface_n_resi': 'trunc_iface_n_resi',
            'binder_iface_plddt': 'trunc_iface_plddt',
            'binder_iface_residues': 'trunc_iface_residues',
            'binder_iface_intracellular': 'trunc_iface_intracellular',
            'binder_resi_start': 'trunc_resi_start',
            'binder_resi_end': 'trunc_resi_end',
            'binder_length': 'trunc_length',
            'binder_total_n_ss': 'trunc_total_n_ss',
            'binder_n_structured_ss': 'trunc_n_structured_ss',
            'binder_structured_percent': 'trunc_structured_percent',
            'binder_iface_n_ss': 'trunc_iface_n_ss',
            'binder_iface_structured': 'trunc_iface_structured',
            'binder_iface_ss': 'trunc_iface_ss',
        }
        binder_metrics_trunc_renamed = {}
        for key, new_key in mapping.items():
            if key in binder_metrics_trunc:
                binder_metrics_trunc_renamed[new_key] = binder_metrics_trunc[key]

        # Recompute geodesic metrics for truncated binder.
        try:
            geo_metrics_trunc = compute_geodesic_length(str(trunc_binder_path))
            # Apply the same scaling factor as before (multiply by 1000)
            if geo_metrics_trunc['norm_geodesic_length'] is not None:
                geo_metrics_trunc['norm_geodesic_length'] = geo_metrics_trunc['norm_geodesic_length'] * 1000
        except Exception as e:
            print(f"Error computing geodesic length for {trunc_binder_path}: {e}")
            geo_metrics_trunc = {'abs_geodesic_length': None, 'norm_geodesic_length': None}
        geo_metrics_trunc_renamed = {}
        if "abs_geodesic_length" in geo_metrics_trunc:
            geo_metrics_trunc_renamed["abs_geodesic_length_trunc"] = geo_metrics_trunc["abs_geodesic_length"]
        if "norm_geodesic_length" in geo_metrics_trunc:
            geo_metrics_trunc_renamed["norm_geodesic_length_trunc"] = geo_metrics_trunc["norm_geodesic_length"]

        # Compute new difference metrics between truncated and original binder interface metrics.
        # Original binder metrics are assumed to be in the CSV row.
        d_iface_n_resi = None
        d_iface_n_ss = None
        d_iface_structured = None
        try:
            orig_iface_n_resi = float(row.get("binder_iface_n_resi", 0))
            trunc_iface_n_resi = float(binder_metrics_trunc_renamed.get("trunc_iface_n_resi", 0))
            d_iface_n_resi = trunc_iface_n_resi - orig_iface_n_resi
        except Exception:
            pass
        try:
            orig_iface_n_ss = float(row.get("binder_iface_n_ss", 0))
            trunc_iface_n_ss = float(binder_metrics_trunc_renamed.get("trunc_iface_n_ss", 0))
            d_iface_n_ss = trunc_iface_n_ss - orig_iface_n_ss
        except Exception:
            pass
        try:
            orig_iface_structured = float(row.get("binder_iface_structured", 0))
            trunc_iface_structured = float(binder_metrics_trunc_renamed.get("trunc_iface_structured", 0))
            d_iface_structured = trunc_iface_structured - orig_iface_structured
        except Exception:
            pass

        # Combine new metrics into one dictionary.
        new_metrics = {}
        new_metrics["clashes_heavy_strictness1_trunc"] = clashes1
        new_metrics["clashes_heavy_strictness0.75_trunc"] = clashes075
        new_metrics.update(sasa_metrics)
        new_metrics.update(binder_metrics_trunc_renamed)
        new_metrics.update(geo_metrics_trunc_renamed)
        new_metrics["d_iface_n_resi"] = d_iface_n_resi
        new_metrics["d_iface_n_ss"] = d_iface_n_ss
        new_metrics["d_iface_structured"] = d_iface_structured

        # ---- NEW: Extract EvoEF2_score from trunc_binder_path filename ----
        try:
            basename = trunc_binder_path.name  # e.g., "Q8NA31-F1-dom-01_A_1465_101_186_-451.pdb"
            evoef2_str = basename.rsplit("_", 1)[1].replace(".pdb", "")
            new_metrics["EvoEF2_score"] = float(evoef2_str)
        except Exception as e:
            new_metrics["EvoEF2_score"] = None

        # ---- NEW: Compute norm_EvoEF2_score as EvoEF2_score divided by trunc_length ----
        try:
            # Ensure both values are available and trunc_length is non-zero.
            trunc_length = new_metrics.get("trunc_length")
            evoef2_score = new_metrics.get("EvoEF2_score")
            if evoef2_score is not None and trunc_length and trunc_length != 0:
                new_metrics["norm_EvoEF2_score"] = evoef2_score / trunc_length
            else:
                new_metrics["norm_EvoEF2_score"] = None
        except Exception:
            new_metrics["norm_EvoEF2_score"] = None
        
        # ---- NEW: Compute norm_trunc_sasa as trunc_sasa divided by trunc_length ----
        try:
            trunc_sasa = new_metrics.get("trunc_sasa")
            # Reuse trunc_length from above
            if trunc_sasa is not None and trunc_length and trunc_length != 0:
                new_metrics["norm_trunc_sasa"] = trunc_sasa / trunc_length
            else:
                new_metrics["norm_trunc_sasa"] = None
        except Exception:
            new_metrics["norm_trunc_sasa"] = None

        # Merge the new metrics with the original row's data
        row_dict = row.to_dict()
        row_dict.update(new_metrics)
        results.append(row_dict)

    return pd.DataFrame(results)

# ------------------ Main ------------------
if __name__ == "__main__":
    parser = ArgumentParser(
        description="Process an input CSV to compute metrics for truncated binder PDB files."
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to the input CSV file (created by postprocess_masif_hits.py)."
    )
    parser.add_argument(
        "-o", "--out_csv_file",
        type=Path,
        required=True,
        help="Path to the output CSV file."
    )
    parser.add_argument(
        "--ligand",
        type=str,
        default=None,
        help="Ligand definition in the format 'CHAIN_RESNAME', e.g., 'A_LIG'."
    )
    parser.add_argument(
        "--sulfur",
        type=str,
        default=None,
        help="Atom name of the sulfur. (Optional)"
    )
    
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print debug information (e.g., deepTMHMM and STRIDE output)."
    )
    args = parser.parse_args()

    # Note: In this script we do not use parse_score_file since scores are already in the input CSV.

    if args.ligand is not None:
        try:
            ligand_chain, ligand_name = args.ligand.split("_", 1)
            ligand = {"chain": ligand_chain, "name": ligand_name}
        except ValueError:
            print("Error: Ligand must be in the format 'CHAIN_RESNAME', e.g., 'A_LIG'.")
            sys.exit(1)
    else:
        ligand = None

    df_trunc = process_truncated_results(
        input_csv=args.input,
        ligand_def=ligand,
        sulfur_name=args.sulfur,
        debug=args.debug
    )
    df_trunc.to_csv(args.out_csv_file, index=False)
    print(f"Results saved to {args.out_csv_file}")
#!/usr/bin/env python3
"""
postprocess_masif_mimicry.py

This script processes a list of PDB files (binder structures) to compute various structural and interface metrics 
for MaSIF grid search results. The metrics computed include clash counts, SASA metrics, interface metrics,
geodesic lengths, and a ligand descriptor distance score computed by an external script.

Usage:
    python3 /scratch/ymeng/MaSIF_surface_mimicry/postprocess_masif_mimicry.py \
        --pdb_paths /scratch/ymeng/MaSIF_surface_mimicry/6H0F_C_B_mimicry_binders/A0A0A0MRZ7-F1-dom-01.pdb \
        --target_pdb_path /scratch/ymeng/MaSIF_surface_mimicry/Preprocess/input/6H0F_B_B_Y70_502/6H0F_B.pdb \
        --query_preprocessed_root /scratch/ymeng/MaSIF_surface_mimicry/Preprocess \
        --out_csv_file /scratch/ymeng/MaSIF_surface_mimicry/out.csv \
        --ligand B_Y70 \
        --debug

Author: Yanxiang Meng (10.3.2025)
"""

from pathlib import Path
from argparse import ArgumentParser
import subprocess
import math  
from rdkit import Chem
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.DSSP import DSSP
from Bio.PDB.Polypeptide import three_to_one  # for sequence extraction
from Bio import pairwise2  # for alignment
from Bio.PDB.NeighborSearch import NeighborSearch  # for efficient spatial searching
from Bio.PDB import StructureBuilder
import yaml
import os

# Load configuration from config/config.yaml
CONFIG_PATH = os.path.join(os.path.dirname(__file__), "config.yaml")
with open(CONFIG_PATH, 'r') as f:
    config = yaml.safe_load(f)

paths_config = config['paths']
deeptmhmm_dir = Path(paths_config['deeptmhmm_dir'])
python_path = paths_config['python_path']
stride_exec = paths_config['stride_exec']
masif_grid_search_scripts_remote = paths_config['masif_grid_search_scripts_remote']

filenames_config = config['filenames']
geosedic_py_path = os.path.join(masif_grid_search_scripts_remote, filenames_config['geosedic_py'])
lig_desc_dist_score_py = os.path.join(masif_grid_search_scripts_remote, filenames_config['lig_desc_dist_score_py'])

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
        ligand_residue = extract_ligand_residue(target_struct, ligand_chain=ligand["chain"], ligand_name=ligand["name"])
        copy_residues([ligand_residue], out_struct, output_chain="L")
    return out_struct

def merge_binder_ligand(binder_struct, ligand_residue):
    """Builds a structure containing only the binder and the ligand."""
    builder = StructureBuilder.StructureBuilder()
    builder.init_structure('binder_ligand_complex')
    builder.init_model(0)
    builder.init_chain("B")  # binder
    builder.init_chain("L")  # ligand
    new_struct = builder.get_structure()[0]
    copy_residues(binder_struct, new_struct, output_chain="B")
    copy_residues([ligand_residue], new_struct, output_chain="L")
    return new_struct

def compute_sasa_values(target_pdb, binder_pdb, ligand_def: dict = None):
    out = {}

    target = PDBParser(QUIET=True).get_structure('', target_pdb)[0]
    binder = PDBParser(QUIET=True).get_structure('', binder_pdb)[0]
    whole_complex = merge_structures(target, binder, ligand_def)

    # Target
    ShrakeRupley().compute(target, level="R")
    out['target_sasa'] = sum(r.sasa for r in target.get_residues())

    # Binder
    ShrakeRupley().compute(binder, level="R")
    out['binder_sasa'] = sum(r.sasa for r in binder.get_residues())

    # Full Complex (target + binder [+ ligand])
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
    if ligand_def is not None:
        ligand_res = extract_ligand_residue(target_struct, ligand_def["chain"], ligand_def["name"])
        contact_atoms = [atom for atom in ligand_res if atom.element != 'H']
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
    """Return binder residues that have any heavy atom within the cutoff distance to any contact atom.
       This version uses NeighborSearch for efficiency."""
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
    """Compute average B-factor for interface residues using only the alpha carbon (CA) atoms."""
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
                            aa = three_to_one(res.get_resname())
                        except Exception:
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




def compute_geodesic_lengths(pdb_path):
    """Compute geodesic length metrics by calling an external script."""
    try:
        result = subprocess.run(
            [python_path, geosedic_py_path, "--input", str(pdb_path), "--abs", "--norm"],
            capture_output=True, text=True, check=True
        )
        lines = result.stdout.strip().splitlines()
        if len(lines) >= 2:
            return {
                'abs_geosedic_length': float(lines[0]),
                'norm_geosedic_length': float(lines[1]),
            }
        else:
            return {'abs_geosedic_length': None, 'norm_geosedic_length': None}
    except Exception as e:
        print(f"Error calling geodesic_length.py for {pdb_path}: {e}")
        return {'abs_geosedic_length': None, 'norm_geosedic_length': None}

# ------------------ Main Processing Function ------------------
def process_results(
    pdb_paths, 
    target_pdb, 
    query_preprocessed_root,
    ligand_def: dict = None, 
    sulfur_name: str = None, 
    recompute_clashes: bool = True, 
    compute_sasa_flag: bool = True,
    debug: bool = False
):
    results = []
    # Cache parsed structures to avoid re-parsing the same PDB file multiple times.
    parsed_structures = {}
    parser = PDBParser(QUIET=True)
    
    # Parse the target structure from the user-provided path.
    target_pdb = Path(target_pdb).resolve()
    try:
        target_struct = parser.get_structure('', str(target_pdb))[0]
    except Exception as e:
        print(f"Error parsing target structure {target_pdb}: {e}")
        return pd.DataFrame(results)
    
    for pdb_path in tqdm(pdb_paths, desc="Processing PDB files"):
        pdb_path = Path(pdb_path).resolve()
        if not pdb_path.is_file():
            print(f"Warning: {pdb_path} is not a valid file. Skipping.")
            continue

        # Parse binder structure.
        if pdb_path not in parsed_structures:
            try:
                parsed_structures[pdb_path] = parser.get_structure('', str(pdb_path))[0]
            except Exception as e:
                print(f"Error parsing binder structure {pdb_path}: {e}")
                continue
        binder_struct = parsed_structures[pdb_path]

        match_info = {}
        match_info['matched_protein_path'] = str(pdb_path)
        # Add target_path and ligand columns to the output.
        match_info['target_path'] = str(target_pdb)
        match_info['ligand'] = f"{ligand_def['chain']}_{ligand_def['name']}" if ligand_def is not None else None

        if recompute_clashes:
            match_info["clashes_heavy_strictness1"] = count_clashes(str(target_pdb), str(pdb_path), strictness=1.0)
            match_info["clashes_heavy_strictness0.75"] = count_clashes(str(target_pdb), str(pdb_path), strictness=0.75)
        if compute_sasa_flag:
            match_info.update(compute_sasa_values(str(target_pdb), str(pdb_path), ligand_def=ligand_def))

        binder_metrics = compute_binder_interface_metrics(binder_struct, target_struct, pdb_path, ligand_def, debug=debug)
        match_info.update(binder_metrics)

        if sulfur_name is not None and ligand_def is not None:
            match_info.update({
                **(lambda: 
                    # Calculate distance from query sulfur atom to closest CB in binder structure.
                    # (Retaining the same logic as before.)
                    (lambda query:
                        (
                            lambda query_struct, target_struct:
                                (
                                    lambda c_betas:
                                        {
                                            "closest_c_beta": f"{c_betas[np.argmin(np.sqrt(np.sum((query['coord'] - np.stack([a.get_coord() for a in c_betas]))**2, axis=1)))]}" if c_betas else None,
                                            "closest_c_beta_distance": np.min(np.sqrt(np.sum((query['coord'] - np.stack([a.get_coord() for a in c_betas]))**2, axis=1))) if c_betas else None
                                        }
                                    )( [atom for atom in target_struct.get_atoms() if atom.name == "CB"] )
                                )( 
                                    PDBParser(QUIET=True).get_structure("", str(target_pdb))[0]
                                )
                        )
                    )({"pdb": str(target_pdb), "chain": ligand_def["chain"], "resname": ligand_def["name"], "atom_name": sulfur_name})
                
            })
            
        # --- Compute geodesic lengths ---
        geodesic_metrics = compute_geodesic_lengths(pdb_path)
        match_info.update(geodesic_metrics)

        results.append(match_info)
    return pd.DataFrame(results)

if __name__ == "__main__":
    from sys import exit
    parser = ArgumentParser(
        description="Process PDB files to compute clashes, SASA values, and interface metrics."
    )
    parser.add_argument(
        "--pdb_paths", 
        nargs='+', 
        type=Path, 
        required=True, 
        help="List of absolute paths to binder .pdb files."
    )
    parser.add_argument(
        "--target_pdb_path", 
        type=Path, 
        required=True, 
        help="Absolute path to the target .pdb file."
    )
    parser.add_argument(
        "--query_preprocessed_root", 
        type=Path, 
        required=True, 
        help="Directory containing the preprocessed query data."
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
        help="Atom name of the sulfur. Needs to be found in `ligand`."
    )
    parser.add_argument(
        "--debug", 
        action="store_true", 
        help="Print debug information (DeepTMHMM and STRIDE output)."
    )
    args = parser.parse_args()

    # Parse ligand if provided.
    if args.ligand is not None:
        try:
            ligand_chain, ligand_name = args.ligand.split("_", 1)
            ligand = {"chain": ligand_chain, "name": ligand_name}
        except ValueError:
            print("Error: Ligand must be in the format 'CHAIN_RESNAME', e.g., 'A_LIG'.")
            exit(1)
    else:
        ligand = None

    df = process_results(
        pdb_paths=args.pdb_paths, 
        target_pdb=args.target_pdb_path,
        query_preprocessed_root=args.query_preprocessed_root,
        ligand_def=ligand, 
        sulfur_name=args.sulfur,
        debug=args.debug
    )
    df.to_csv(args.out_csv_file, index=False)
    print(f"Results saved to {args.out_csv_file}")
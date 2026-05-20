"""Compute post-processing metrics for deduplicated mimicry search hits."""

import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBIO
from tqdm import tqdm

from masif_mimicry.postprocess.metrics.clashes import count_clashes
from masif_mimicry.postprocess.metrics.interface import compute_binder_interface_metrics
from masif_mimicry.postprocess.metrics.sasa import compute_sasa_values
from masif_mimicry.postprocess.structures import maybe_load_structure
from masif_mimicry.utils.transforms import get_transformed_struct_from_row


def parse_ligand_def(ligand_def: str) -> dict:
    """Parse --ligand CHAIN_RESNAME (e.g. A_021 -> chain A, hetero resname 021)."""
    if "_" not in ligand_def:
        raise ValueError(f"--ligand must be CHAIN_RESNAME, got {ligand_def!r}")
    chain, resname = ligand_def.split("_", 1)
    return {"chain": chain, "name": resname}


def load_patch_coord(target_preprocess_dir, target_name, center_vix):
    precomputation_dir = Path(
        target_preprocess_dir,
        "data_preparation",
        "04b-precomputation_12A",
        "precomputation",
        target_name,
    )
    vertex_indices = np.load(Path(precomputation_dir, "p1_list_indices.npy"), allow_pickle=True)
    vertex_coord = np.stack(
        [np.load(Path(precomputation_dir, f"p1_{dim}.npy")) for dim in ["X", "Y", "Z"]],
        axis=1,
    )
    return vertex_coord[vertex_indices[center_vix]]


def structure_to_temp_pdb(struct):
    fd, path = tempfile.mkstemp(suffix=".pdb")
    os.close(fd)
    io = PDBIO()
    io.set_structure(struct)
    io.save(path)
    return path


SASA_PLACEHOLDER = {
    "target_unbound_sasa": None,
    "binder_unbound_sasa": None,
    "ligand_unbound_sasa": None,
    "target_ligand_sasa": None,
    "binder_ligand_sasa": None,
    "complex_sasa": None,
    "target_in_complex_sasa": None,
    "binder_in_complex_sasa": None,
    "ligand_in_complex_sasa": None,
    "ligand_in_lb_sasa": None,
    "binder_in_lb_sasa": None,
    "ligand_in_tl_sasa": None,
    "target_in_tl_sasa": None,
    "target_buried_in_complex": None,
    "binder_buried_in_complex": None,
    "ligand_buried_in_complex": None,
    "target_buried_in_tb": None,
    "binder_buried_in_tb": None,
    "ligand_buried_in_lb": None,
    "binder_buried_in_lb": None,
    "ligand_buried_in_tl": None,
    "target_buried_in_tl": None,
    "ligand_iface_contribution": None,
    "delta_sasa": None,
}


def process_results_mimicry(
    df,
    target_pdb,
    database_dir,
    target_preprocess_dir,
    ligand_def,
    out_csv_file=None,
):
    """
    Loop over deduplicated rows and append metric columns to each row.
    Preserves existing CSV columns (including flattened_transform).
    """
    target_pdb = Path(target_pdb)
    ligand = parse_ligand_def(ligand_def)
    results = []
    first_write = True

    print(f"Postprocessing {len(df)} rows...")
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="postprocess"):
        match_info = row.to_dict()
        p1_id = match_info.get("P1_id", idx)
        binder_tmp = None

        try:
            matched_struct = get_transformed_struct_from_row(row, database_dir)
        except Exception as e:
            print(f"[{p1_id}] Error building transformed binder: {e}", flush=True)
            matched_struct = None

        try:
            target_struct = maybe_load_structure(target_pdb)
        except Exception as e:
            print(f"[{p1_id}] Error loading target structure: {e}", flush=True)
            target_struct = None

        if matched_struct is not None:
            try:
                binder_tmp = structure_to_temp_pdb(matched_struct)
                match_info["clashes_heavy_strictness1"] = count_clashes(
                    target_pdb, binder_tmp, strictness=1.0
                )
                match_info["clashes_heavy_strictness0.75"] = count_clashes(
                    target_pdb, binder_tmp, strictness=0.75
                )
            except Exception as e:
                print(f"[{p1_id}] Error computing clashes: {e}", flush=True)
                match_info["clashes_heavy_strictness1"] = None
                match_info["clashes_heavy_strictness0.75"] = None
        else:
            match_info["clashes_heavy_strictness1"] = None
            match_info["clashes_heavy_strictness0.75"] = None

        try:
            patch_coord = load_patch_coord(
                target_preprocess_dir,
                match_info["P2_id"],
                int(match_info["P2_source_site"]),
            )
            binder_metrics = compute_binder_interface_metrics(
                matched_struct,
                target_struct,
                target_patch_coord=patch_coord,
            )
            match_info["target_iface_resi"] = binder_metrics.get("target_iface_resi")
            match_info["matched_iface_resi"] = binder_metrics.get("matched_iface_resi")
            match_info["matched_iface_n_resi"] = binder_metrics.get("matched_iface_n_resi")
            match_info["matched_iface_plddt"] = binder_metrics.get("matched_iface_plddt")
        except Exception as e:
            print(f"[{p1_id}] Error computing binder interface metrics: {e}", flush=True)
            match_info["target_iface_resi"] = None
            match_info["matched_iface_resi"] = None
            match_info["matched_iface_n_resi"] = None
            match_info["matched_iface_plddt"] = None

        try:
            if binder_tmp is None and matched_struct is not None:
                binder_tmp = structure_to_temp_pdb(matched_struct)
            sasa_metrics = compute_sasa_values(target_pdb, binder_tmp, ligand_def=ligand)
            match_info.update(sasa_metrics)
        except Exception as e:
            print(f"[{p1_id}] Error computing SASA values: {e}", flush=True)
            match_info.update(SASA_PLACEHOLDER)
        finally:
            if binder_tmp and os.path.isfile(binder_tmp):
                os.remove(binder_tmp)

        if out_csv_file is not None:
            df_row = pd.DataFrame([match_info])
            mode = "w" if first_write else "a"
            df_row.to_csv(out_csv_file, mode=mode, header=first_write, index=False)
            first_write = False
        else:
            results.append(match_info)

    if out_csv_file is not None:
        return None
    return pd.DataFrame(results)

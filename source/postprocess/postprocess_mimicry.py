import os
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBIO
from tqdm import tqdm

if __name__ == "__main__":
    basedir = Path(__file__).resolve().parent.parent
    sys.path.append(str(basedir))

from utils import get_transformed_struct_from_row, resolve_database_paths
from postprocess.postprocess_utils import maybe_load_structure
from postprocess.metrics.clashes import count_clashes
from postprocess.metrics.sasa import compute_sasa_values
from postprocess.metrics.interface import compute_binder_interface_metrics


"""
Process mimicry search outputs: discover per-seed CSVs under --target_run_dir,
deduplicate to one row per P1, and compute post-processing metrics.

Command-line usage:
    python process_search_outputs_mimicry.py \\
        --target_run_dir /path/to/021structure_C_AB_ \\
        --target_pdb /path/to/021structure_C.pdb \\
        --database_dir /path/to/TED_domainome/output \\
        --target_preprocess_dir /path/to/NUP98 \\
        --ligand A_021 \\
        -o output.csv \\
        [--subset subset_ids.txt]
"""


def parse_ligand_def(ligand_def: str) -> dict:
    """
    Parse --ligand CHAIN_RESNAME (e.g. A_021 -> chain A, hetero resname 021).
    """
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


def _cluster_summary(df):
    work = df.copy()
    work["cluster_size"] = pd.to_numeric(work["cluster_size"], errors="coerce").fillna(0)
    work["cluster_mean_rmsd"] = pd.to_numeric(work["cluster_mean_rmsd"], errors="coerce")
    summary = (
        work.groupby("cluster_id", dropna=False)
        .agg(
            cluster_size=("cluster_size", "max"),
            cluster_mean_rmsd=("cluster_mean_rmsd", "min"),
        )
        .reset_index()
    )
    return summary


def select_representative_row(df):
    """Pick one row per seed CSV using cluster_size / cluster_mean_rmsd / MaSIF-score."""
    summary = _cluster_summary(df)
    max_size = summary["cluster_size"].max()
    top_clusters = summary[summary["cluster_size"] == max_size]

    if len(top_clusters) == 1:
        chosen_cluster = top_clusters["cluster_id"].iloc[0]
    else:
        rmsd = top_clusters["cluster_mean_rmsd"].fillna(np.inf)
        min_rmsd = rmsd.min()
        chosen_cluster = top_clusters.loc[rmsd == min_rmsd, "cluster_id"].iloc[0]

    candidates = df[df["cluster_id"] == chosen_cluster]
    masif = pd.to_numeric(candidates["MaSIF-score"], errors="coerce")
    best_score = masif.max()
    return candidates.loc[masif == best_score].iloc[0]


def discover_deduplicated_rows(target_run_dir):
    target_run_dir = Path(target_run_dir)
    rows = []
    skip_names = {"target_vert"}

    for p1_dir in sorted(target_run_dir.iterdir()):
        if not p1_dir.is_dir() or p1_dir.name in skip_names:
            continue
        if p1_dir.name.endswith(".pdb") or p1_dir.name.endswith(".json"):
            continue

        p1 = p1_dir.name
        csvs = sorted(p1_dir.glob(f"{p1}_p1_to_*.csv"))
        if len(csvs) == 0:
            print(f"No hit found for {p1}", flush=True)
            continue
        if len(csvs) > 1:
            warnings.warn(
                f"Skipping {p1}: expected one *_p1_to_*.csv, found {len(csvs)}: "
                + ", ".join(c.name for c in csvs)
            )
            continue

        hit_df = pd.read_csv(csvs[0])
        if len(hit_df) == 0:
            print(f"No hit found for {p1} (empty CSV)", flush=True)
            continue
        rows.append(select_representative_row(hit_df))

    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).reset_index(drop=True)


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


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Post-process mimicry search CSV outputs")
    parser.add_argument(
        "--target_run_dir",
        type=Path,
        required=True,
        help="Search run directory containing per-P1 subfolders and CSVs",
    )
    parser.add_argument(
        "--target_pdb",
        type=Path,
        required=True,
        help="Path to target chain-C PDB (e.g. 021structure_C.pdb)",
    )
    parser.add_argument(
        "--database_dir",
        type=Path,
        required=True,
        help="MaSIF database root (e.g. TED_domainome/output); PDBs under data_preparation/01-benchmark_pdbs/",
    )
    parser.add_argument(
        "--target_preprocess_dir",
        type=Path,
        required=True,
        help="MaSIF preprocess directory for the target (e.g. data/NUP98)",
    )
    parser.add_argument("-o", "--out_csv_file", type=Path, required=True, help="Combined output CSV")
    parser.add_argument("--subset", type=Path, default=None, help="Optional P1_id list (one per line)")
    parser.add_argument("--ligand", type=str, required=True, help="Ligand name and chain, e.g. 'A_021'")

    args = parser.parse_args()

    database_dir = os.path.abspath(os.path.expanduser(args.database_dir))
    _, db_prep = resolve_database_paths(database_dir)
    if not os.path.isdir(db_prep):
        print(f"Error: data_preparation directory not found: {db_prep}", file=sys.stderr)
        sys.exit(1)

    df = discover_deduplicated_rows(args.target_run_dir)
    if df.empty:
        print("No hits to post-process.", flush=True)
        sys.exit(0)

    if args.subset is not None:
        with open(args.subset) as f:
            subset_list = {line.strip() for line in f if line.strip()}
        df = df[df["P1_id"].isin(subset_list)]

    process_results_mimicry(
        df,
        target_pdb=args.target_pdb,
        database_dir=database_dir,
        target_preprocess_dir=args.target_preprocess_dir,
        ligand_def=args.ligand,
        out_csv_file=args.out_csv_file,
    )
    print(f"Results written to {args.out_csv_file}")

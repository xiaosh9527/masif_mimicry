"""Discover and deduplicate per-seed mimicry search CSV outputs."""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd

SKIP_DIR_NAMES = frozenset({"target_vert"})


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
    """Load one representative row per P1 from *_p1_to_*.csv files under target_run_dir."""
    target_run_dir = Path(target_run_dir)
    rows = []

    for p1_dir in sorted(target_run_dir.iterdir()):
        if not p1_dir.is_dir() or p1_dir.name in SKIP_DIR_NAMES:
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

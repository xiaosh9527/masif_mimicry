"""Discover and deduplicate per-seed mimicry search CSV outputs."""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from masif_mimicry.postprocess.metrics.patch_consistency import (
    load_seed_surface_coords,
    load_target_surface_coords,
    parse_transform,
    patch_consistency_score,
)

SKIP_DIR_NAMES = frozenset({"target_vert"})
CLUSTER_METHODS = ("max_score", "medoid")


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


def _choose_cluster_id(df):
    """Select cluster_id by largest cluster_size, then lowest cluster_mean_rmsd."""
    summary = _cluster_summary(df)
    max_size = summary["cluster_size"].max()
    top_clusters = summary[summary["cluster_size"] == max_size]

    if len(top_clusters) == 1:
        return top_clusters["cluster_id"].iloc[0]

    rmsd = top_clusters["cluster_mean_rmsd"].fillna(np.inf)
    min_rmsd = rmsd.min()
    return top_clusters.loc[rmsd == min_rmsd, "cluster_id"].iloc[0]


def n_p2_source_sites_in_cluster(df, cluster_id):
    """Count unique P2_source_site values among hits in cluster_id."""
    cluster_hits = df[df["cluster_id"] == cluster_id]
    if len(cluster_hits) == 0 or "P2_source_site" not in cluster_hits.columns:
        return 0
    return int(cluster_hits["P2_source_site"].nunique())


def _representative_index(candidates, cluster_method):
    """Index of the cluster representative row per ``cluster_method``.

    ``max_score``: highest MaSIF-score. ``medoid``: most geometrically central
    pose, i.e. lowest ``cluster_mean_rmsd`` (mean backbone RMSD to the other
    cluster members, precomputed in :mod:`masif_mimicry.utils.clustering`). The
    medoid falls back to max_score when ``cluster_mean_rmsd`` is unavailable.
    """
    if cluster_method == "medoid":
        rmsd = pd.to_numeric(candidates["cluster_mean_rmsd"], errors="coerce")
        if rmsd.notna().any():
            return rmsd.idxmin()
    masif = pd.to_numeric(candidates["MaSIF-score"], errors="coerce")
    return masif.idxmax()


def search_hit_csv_path(
    target_run_dir, p1_id, p2_id, seed_ppi_id="p1", target_ppi_id="p1"
):
    """Predicted search-hit CSV path from mimicry search output layout."""
    target_run_dir = Path(target_run_dir)
    return target_run_dir / p1_id / f"{p1_id}_{seed_ppi_id}_to_{p2_id}_{target_ppi_id}.csv"


def select_representative_row(
    df, seed_surf_coord=None, target_surf_coord=None, cluster_method="max_score"
):
    """Pick one row per seed CSV from the chosen cluster.

    ``cluster_method`` selects the representative within the cluster: ``max_score``
    (highest MaSIF-score) or ``medoid`` (most geometrically central pose). The same
    representative pose is used downstream and as the reference for the
    patch-consistency metrics.

    When ``seed_surf_coord`` and ``target_surf_coord`` are provided, also attach
    patch-consistency metrics measuring how well the representative pose explains
    the other patch matches in the cluster.
    """
    if cluster_method not in CLUSTER_METHODS:
        raise ValueError(f"cluster_method must be one of {CLUSTER_METHODS}, got {cluster_method!r}")

    chosen_cluster = _choose_cluster_id(df)
    candidates = df[df["cluster_id"] == chosen_cluster]
    n_p2 = n_p2_source_sites_in_cluster(df, chosen_cluster)

    row = candidates.loc[_representative_index(candidates, cluster_method)].copy()
    row["n_P2_source_site"] = n_p2

    try:
        reference_transform = parse_transform(row["flattened_transform"])
    except Exception:
        reference_transform = None
    consistency = patch_consistency_score(
        candidates, seed_surf_coord, target_surf_coord, reference_transform=reference_transform
    )
    for key, value in consistency.items():
        row[key] = value
    return row


def _coord_loaders(database_dir, target_preprocess_dir):
    target_coord_cache = {}

    def _target_coords(target_name):
        if target_preprocess_dir is None:
            return None
        if target_name not in target_coord_cache:
            try:
                target_coord_cache[target_name] = load_target_surface_coords(
                    target_preprocess_dir, target_name
                )
            except Exception as e:
                warnings.warn(f"Could not load target surface coords for {target_name}: {e}")
                target_coord_cache[target_name] = None
        return target_coord_cache[target_name]

    def _seed_coords(seed_name):
        if database_dir is None:
            return None
        try:
            return load_seed_surface_coords(database_dir, seed_name)
        except Exception as e:
            warnings.warn(f"Could not load seed surface coords for {seed_name}: {e}")
            return None

    return _seed_coords, _target_coords


def _deduplicated_row_for_hit_csv(
    p1,
    csv_path,
    *,
    database_dir=None,
    target_preprocess_dir=None,
    cluster_method="max_score",
    seed_coords=None,
    target_coords=None,
):
    hit_df = pd.read_csv(csv_path)
    if len(hit_df) == 0:
        print(f"No hit found for {p1} (empty CSV)", flush=True)
        return None

    if seed_coords is None or target_coords is None:
        seed_coords, target_coords = _coord_loaders(database_dir, target_preprocess_dir)
    seed_surf_coord = seed_coords(str(hit_df["P1_id"].iloc[0]))
    target_surf_coord = target_coords(str(hit_df["P2_id"].iloc[0]))
    return select_representative_row(
        hit_df, seed_surf_coord, target_surf_coord, cluster_method=cluster_method
    )


def discover_deduplicated_rows_for_subset(
    target_run_dir,
    p1_ids,
    p2_id,
    seed_ppi_id="p1",
    target_ppi_id="p1",
    database_dir=None,
    target_preprocess_dir=None,
    cluster_method="max_score",
):
    """Load one representative row per P1_id listed in a subset file.

    Uses the predictable search CSV path
    ``{target_run_dir}/{P1_id}/{P1_id}_{seed_ppi_id}_to_{p2_id}_{target_ppi_id}.csv``.
    """
    target_run_dir = Path(target_run_dir)
    seed_coords, target_coords = _coord_loaders(database_dir, target_preprocess_dir)
    rows = []

    for p1 in p1_ids:
        p1 = str(p1).strip()
        if not p1:
            continue

        csv_path = search_hit_csv_path(
            target_run_dir, p1, p2_id, seed_ppi_id=seed_ppi_id, target_ppi_id=target_ppi_id
        )
        if not csv_path.is_file():
            print(f"No hit found for {p1}", flush=True)
            continue

        row = _deduplicated_row_for_hit_csv(
            p1,
            csv_path,
            database_dir=database_dir,
            target_preprocess_dir=target_preprocess_dir,
            cluster_method=cluster_method,
            seed_coords=seed_coords,
            target_coords=target_coords,
        )
        if row is not None:
            rows.append(row)

    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).reset_index(drop=True)


def discover_deduplicated_rows(
    target_run_dir, database_dir=None, target_preprocess_dir=None, cluster_method="max_score"
):
    """Load one representative row per P1 from *_p1_to_*.csv files under target_run_dir.

    ``cluster_method`` (``max_score`` or ``medoid``) selects the representative row
    within each chosen cluster during deduplication.

    When ``database_dir`` and ``target_preprocess_dir`` are given, surface vertex
    coordinates are loaded so each representative row is annotated with
    patch-consistency metrics (see :mod:`masif_mimicry.postprocess.metrics.patch_consistency`).
    """
    target_run_dir = Path(target_run_dir)
    rows = []
    seed_coords, target_coords = _coord_loaders(database_dir, target_preprocess_dir)

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

        row = _deduplicated_row_for_hit_csv(
            p1,
            csvs[0],
            database_dir=database_dir,
            target_preprocess_dir=target_preprocess_dir,
            cluster_method=cluster_method,
            seed_coords=seed_coords,
            target_coords=target_coords,
        )
        if row is not None:
            rows.append(row)

    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).reset_index(drop=True)

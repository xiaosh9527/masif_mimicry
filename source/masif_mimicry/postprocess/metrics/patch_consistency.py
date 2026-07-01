"""Patch-consistency score: agreement of a cluster's hits with its representative pose.

Within a structural cluster, the representative hit (highest MaSIF-score) defines a
rigid transform ``T_r`` that maps the seed surface into the target frame. A genuine
binding mode implies that *every* member's own patch correspondence
(``P1_source_site`` -> ``P2_source_site``) is also satisfied by ``T_r``: transforming
the seed's patch vertex by ``T_r`` should land near the target's patch vertex.

We reward (a) covering more distinct target regions and (b) tighter agreement, by
deduplicating on ``P2_source_site`` and summing a smooth distance kernel (TM-score
style) over the distinct target regions.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from masif_mimicry.config.paths import resolve_database_paths

_PRECOMP_SUBDIR = ("04b-precomputation_12A", "precomputation")


def _load_surface_coords(precomputation_root, name):
    """Stack per-vertex X/Y/Z arrays into an (n_vertices, 3) array."""
    coords = np.stack(
        [np.load(Path(precomputation_root, name, f"p1_{dim}.npy")) for dim in ("X", "Y", "Z")],
        axis=1,
    )
    return coords


def load_seed_surface_coords(database_dir, seed_name):
    """Seed (P1) surface vertex coordinates from the MaSIF database precomputation."""
    _, db_prep = resolve_database_paths(str(database_dir))
    return _load_surface_coords(Path(db_prep, *_PRECOMP_SUBDIR), seed_name)


def load_target_surface_coords(target_preprocess_dir, target_name):
    """Target (P2) surface vertex coordinates from the target preprocess tree."""
    root = Path(target_preprocess_dir, "data_preparation", *_PRECOMP_SUBDIR)
    return _load_surface_coords(root, target_name)


def parse_transform(flattened_transform):
    """Parse a comma-separated flattened 4x4 transform into a (4, 4) array."""
    return np.array(list(map(float, str(flattened_transform).split(",")))).reshape(4, 4)


PATCH_CONSISTENCY_COLUMNS = (
    "patch_consistency_score",
    "patch_consistency_mean",
    "n_consistent_P2_site",
    "patch_consistency_mean_dist",
)


def _empty_result():
    return {
        "patch_consistency_score": np.nan,
        "patch_consistency_mean": np.nan,
        "n_consistent_P2_site": np.nan,
        "patch_consistency_mean_dist": np.nan,
    }


def patch_consistency_score(
    cluster_df,
    seed_surf_coord,
    target_surf_coord,
    reference_transform=None,
    d0=5.0,
    inlier_thresh=5.0,
):
    """Quantify how well the representative pose explains a cluster's patch matches.

    Parameters
    ----------
    cluster_df : DataFrame
        Rows of a single chosen cluster. Needs columns ``P1_source_site``,
        ``P2_source_site``, ``MaSIF-score``, ``flattened_transform``.
    seed_surf_coord, target_surf_coord : (n, 3) arrays
        Surface vertex coordinates for the seed (P1) and target (P2), indexed by
        ``P1_source_site`` / ``P2_source_site`` respectively.
    reference_transform : (4, 4) array or None
        Rigid transform of the representative pose (seed -> target frame). When
        None, it is derived from the cluster's highest-MaSIF-score row.
    d0 : float
        Distance scale (A) of the soft-count kernel ``1 / (1 + (d / d0)^2)``.
    inlier_thresh : float
        Distance cutoff (A) for the hard ``n_consistent_P2_site`` count.

    Returns
    -------
    dict with keys in :data:`PATCH_CONSISTENCY_COLUMNS`. All-NaN if the score
    cannot be computed (missing coords, unparseable transform, no valid regions).
    """
    if seed_surf_coord is None or target_surf_coord is None or len(cluster_df) == 0:
        return _empty_result()

    if reference_transform is None:
        masif = pd.to_numeric(cluster_df["MaSIF-score"], errors="coerce")
        if masif.notna().sum() == 0:
            return _empty_result()
        try:
            reference_transform = parse_transform(
                cluster_df.loc[masif.idxmax(), "flattened_transform"]
            )
        except Exception:
            return _empty_result()

    transform = np.asarray(reference_transform, dtype=float)
    R, t = transform[:3, :3], transform[:3, 3]

    p1_sites = pd.to_numeric(cluster_df["P1_source_site"], errors="coerce")
    p2_sites = pd.to_numeric(cluster_df["P2_source_site"], errors="coerce")

    n_seed = len(seed_surf_coord)
    n_target = len(target_surf_coord)

    # Minimum transformed-vs-target distance per distinct target region (P2_source_site).
    best_dist_per_region = {}
    for p1_site, p2_site in zip(p1_sites, p2_sites):
        if np.isnan(p1_site) or np.isnan(p2_site):
            continue
        p1_site, p2_site = int(p1_site), int(p2_site)
        if not (0 <= p1_site < n_seed) or not (0 <= p2_site < n_target):
            continue
        v1_t = seed_surf_coord[p1_site] @ R.T + t
        d = float(np.linalg.norm(v1_t - target_surf_coord[p2_site]))
        prev = best_dist_per_region.get(p2_site)
        if prev is None or d < prev:
            best_dist_per_region[p2_site] = d

    if not best_dist_per_region:
        return _empty_result()

    dists = np.array(list(best_dist_per_region.values()))
    weights = 1.0 / (1.0 + (dists / d0) ** 2)
    score = float(weights.sum())
    return {
        "patch_consistency_score": score,
        "patch_consistency_mean": score / len(dists),
        "n_consistent_P2_site": int((dists <= inlier_thresh).sum()),
        "patch_consistency_mean_dist": float(dists.mean()),
    }

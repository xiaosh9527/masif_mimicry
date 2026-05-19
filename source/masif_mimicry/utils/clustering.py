"""Structural clustering of mimicry search hits by binding mode."""

import numpy as np
from sklearn.cluster import DBSCAN

from masif_mimicry.utils.transforms import get_transformed_struct_from_row


def _backbone_coords(struct, atoms_to_keep=("N", "CA", "C")):
    coords = []
    for atom in struct.get_atoms():
        name = atom.get_name().strip()
        if name in atoms_to_keep:
            coords.append(atom.get_coord())
    if len(coords) == 0:
        raise ValueError("No backbone atoms found in structure")
    return np.stack(coords)


def _pairwise_rmsd(coords1, coords2):
    if coords1.shape != coords2.shape:
        raise ValueError(f"Coordinate shape mismatch: {coords1.shape} vs {coords2.shape}")
    diff = coords1 - coords2
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=-1))))


def structural_clusters(df, rmsd_thresh=5.0, database_dir=None, database_root=None):
    """
    Cluster hits in a mimicry results table by binding mode (RMSD in target frame).

    Expects columns P1_id, flattened_transform.
    Adds cluster_id, cluster_size, cluster_mean_rmsd.

    Rows that cannot be transformed or have incompatible backbone size receive
    cluster_id=-1 and NA for cluster_size / cluster_mean_rmsd.
    """
    if database_dir is None:
        database_dir = database_root
    if database_dir is None:
        raise ValueError("database_dir is required for structural_clusters")

    df = df.copy()
    df["cluster_id"] = -1
    df["cluster_size"] = np.nan
    df["cluster_mean_rmsd"] = np.nan

    if len(df) == 0:
        return df

    positions = []
    valid_indices = []
    for idx, row in df.iterrows():
        try:
            struct = get_transformed_struct_from_row(row, database_dir)
            positions.append(_backbone_coords(struct))
            valid_indices.append(idx)
        except Exception:
            continue

    if len(valid_indices) == 0:
        return df

    if len(valid_indices) == 1:
        idx = valid_indices[0]
        df.loc[idx, "cluster_id"] = 0
        df.loc[idx, "cluster_size"] = 1
        df.loc[idx, "cluster_mean_rmsd"] = 0.0
        return df

    ref_shape = positions[0].shape
    filtered_positions = []
    filtered_indices = []
    for pos, idx in zip(positions, valid_indices):
        if pos.shape == ref_shape:
            filtered_positions.append(pos)
            filtered_indices.append(idx)

    if len(filtered_indices) == 0:
        return df

    if len(filtered_indices) == 1:
        idx = filtered_indices[0]
        df.loc[idx, "cluster_id"] = 0
        df.loc[idx, "cluster_size"] = 1
        df.loc[idx, "cluster_mean_rmsd"] = 0.0
        return df

    n = len(filtered_positions)
    rmsd_vals = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            try:
                rmsd_ij = _pairwise_rmsd(filtered_positions[i], filtered_positions[j])
            except ValueError:
                rmsd_ij = np.inf
            rmsd_vals[i, j] = rmsd_ij
            rmsd_vals[j, i] = rmsd_ij

    labels = DBSCAN(eps=rmsd_thresh, min_samples=2, metric="precomputed").fit_predict(rmsd_vals)
    n_outliers = (labels == -1).sum()
    if n_outliers > 0:
        labels[labels == -1] = np.arange(n_outliers) + labels.max() + 1

    current_cluster_label = 0
    for lb in sorted(set(labels)):
        member_mask = labels == lb
        member_pos = np.where(member_mask)[0]
        cluster_size = int(member_mask.sum())
        sub_rmsd = rmsd_vals[member_mask][:, member_mask]
        mean_rmsd_per_member = sub_rmsd.mean(axis=1)
        for local_i, global_i in enumerate(member_pos):
            idx = filtered_indices[global_i]
            df.loc[idx, "cluster_id"] = current_cluster_label
            df.loc[idx, "cluster_size"] = cluster_size
            df.loc[idx, "cluster_mean_rmsd"] = mean_rmsd_per_member[local_i]
        current_cluster_label += 1

    return df

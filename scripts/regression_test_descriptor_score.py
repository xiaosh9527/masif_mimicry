#!/usr/bin/env python3
"""
One-time regression: patch-only descriptor scoring vs full point-cloud transform.

Compares compute_descriptor_score(transform_patch_coords(...)) against
compute_score_and_clashes after copy.deepcopy(pcd).transform(T) for the same
multidock alignments. Exits 0 if all scores match within tolerance.
"""
import argparse
import copy
import os
import sys

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MASIF_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
sys.path.insert(0, os.path.join(MASIF_ROOT, "masif", "source"))
sys.path.insert(0, os.path.join(MASIF_ROOT, "masif_mimicry", "source"))

from masif_mimicry.config.paths import set_params
from masif_mimicry.search.docking import multidock, transform_patch_coords
from masif_mimicry.search.features import get_features
from masif_mimicry.search.scoring import compute_descriptor_score, compute_score_and_clashes


def score_legacy(P1_all_feats, P1_center, P2_all_feats, P2_center, transformation):
    P1_pcd = copy.deepcopy(P1_all_feats['pcd'])
    P1_pcd.transform(transformation)
    output, _ = compute_score_and_clashes(
        P1_pdb='',
        P1_pcd=P1_pcd,
        P1_descs=P1_all_feats['desc'],
        P1_site=P1_center,
        P1_indices=P1_all_feats['indices'],
        P2_pdb=P2_all_feats['pdb'],
        P2_pcd=P2_all_feats['pcd'],
        P2_descs=P2_all_feats['desc'],
        P2_site=P2_center,
        P2_indices=P2_all_feats['indices'],
        compute_clashes=False,
    )
    return output[1]


def score_deferred(P1_all_feats, P1_center, P2_all_feats, P2_center, transformation):
    P1_patch_coords = transform_patch_coords(
        P1_all_feats['pcd'],
        P1_all_feats['indices'],
        P1_center,
        transformation,
    )
    return compute_descriptor_score(
        P1_patch_coords,
        P1_all_feats['desc'],
        P1_all_feats['indices'],
        P1_center,
        P2_all_feats['pcd'],
        P2_all_feats['desc'],
        P2_all_feats['indices'],
        P2_center,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--database_dir', required=True)
    p.add_argument('--target_preprocess_dir', required=True)
    p.add_argument('--seed_pdb', required=True)
    p.add_argument('--target_pdb', required=True)
    p.add_argument('--target_ppi_id', default='p1', choices=['p1', 'p2'])
    p.add_argument('--seed_ppi_id', default='p1', choices=['p1', 'p2'])
    p.add_argument('--target_site', type=int, default=None, help='P2 center vertex index')
    p.add_argument('--tolerance', type=float, default=1e-12)
    args = p.parse_args()

    params = set_params(
        database_dir=os.path.abspath(args.database_dir),
        target_preprocess_dir=os.path.abspath(args.target_preprocess_dir),
    )
    P1_all_feats = get_features(params, args.seed_pdb, args.seed_ppi_id, source=True, flip_desc=False)
    P2_all_feats = get_features(params, args.target_pdb, args.target_ppi_id, source=False, flip_desc=False)

    P2_center = args.target_site
    if P2_center is None:
        P2_center = int(np.argmax(P2_all_feats['iface'][0]))

    source_pts = np.array([0, len(P1_all_feats['desc']) // 2], dtype=int)
    all_results, _, _, _ = multidock(
        source_pt=source_pts,
        source_pcd=P1_all_feats['pcd'],
        source_patch_idxs=P1_all_feats['indices'],
        source_descs=P1_all_feats['desc'],
        target_pt=P2_center,
        target_pcd=P2_all_feats['pcd'],
        target_patch_idxs=P2_all_feats['indices'],
        target_descs=P2_all_feats['desc'],
        binder_align=False,
    )

    mismatches = []
    for result, P1_center in zip(all_results, source_pts):
        legacy = score_legacy(P1_all_feats, int(P1_center), P2_all_feats, P2_center, result.transformation)
        deferred = score_deferred(P1_all_feats, int(P1_center), P2_all_feats, P2_center, result.transformation)
        diff = abs(legacy - deferred)
        print(f'site P1={P1_center} P2={P2_center}: legacy={legacy:.16f} deferred={deferred:.16f} diff={diff:.3e}')
        if diff > args.tolerance:
            mismatches.append((P1_center, legacy, deferred, diff))

    if mismatches:
        print(f'FAIL: {len(mismatches)} mismatch(es) (tolerance={args.tolerance})', file=sys.stderr)
        sys.exit(1)
    print(f'PASS: {len(all_results)} alignment(s) within tolerance {args.tolerance}')


if __name__ == '__main__':
    main()

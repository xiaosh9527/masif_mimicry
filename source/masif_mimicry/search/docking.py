import numpy as np
from geometry.open3d_import import (
    CorrespondenceCheckerBasedOnDistance,
    CorrespondenceCheckerBasedOnEdgeLength,
    CorrespondenceCheckerBasedOnNormal,
    TransformationEstimationPointToPlane,
    TransformationEstimationPointToPoint,
    registration_icp,
    registration_ransac_based_on_feature_matching,
    RANSACConvergenceCriteria,
)

from masif_mimicry.search.features import ICPConvergenceCriteria, get_patch_geo
from masif_mimicry.utils.transforms import apply_transform


def select_patches(
    P_all_feats: dict,
    downsample_rate: int = 5,
    iface_cutoff: float = 0.0,
    min_patch_num: int = 50,
    top_iface_percent: float = 0.0,
    interface_only: bool = False,
    verbose: bool = True,
) -> tuple:
    """Select patch sites from a MaSIF feature dict."""
    subpcd_coverage = set()
    selected_points_idx = []
    patch_descs = []
    patch_iface = P_all_feats["iface"][0]
    P_indices = P_all_feats["indices"]
    P_interface_points = np.where(P_all_feats["ilabel"] == 1)[0]
    for ii in range(len(P_all_feats["desc"])):
        desc = P_all_feats["desc"][ii]
        if interface_only:
            if ii not in subpcd_coverage and ii in P_interface_points:
                subpcd_coverage.update(P_indices[ii][:downsample_rate])
                selected_points_idx.append(ii)
                patch_descs.append(desc)
        else:
            if ii not in subpcd_coverage and patch_iface[ii] >= iface_cutoff:
                subpcd_coverage.update(P_indices[ii][:downsample_rate])
                selected_points_idx.append(ii)
                patch_descs.append(desc)

    patch_descs = np.array(patch_descs)
    selected_points_idx = np.array(selected_points_idx)

    if len(selected_points_idx) == 0:
        print("No points selected. Please check the parameters.")
        return [], [], []

    if top_iface_percent > 0.0:
        top_iface_num = round(top_iface_percent * len(P_all_feats["desc"]))
        top_iface_num = max(top_iface_num, min_patch_num)
        top_iface_idx = np.argsort(patch_iface[selected_points_idx])[::-1][:top_iface_num]
    else:
        top_iface_idx = np.arange(len(selected_points_idx))

    selected_points_idx = selected_points_idx[top_iface_idx]
    patch_descs = patch_descs[top_iface_idx]

    if verbose:
        if interface_only:
            print(
                f"WARNING: Exhausitive alignment mode. This will go through "
                f"{len(selected_points_idx)} interface points and will take some time..."
            )
        else:
            print(
                f"Selected {len(selected_points_idx)} points with interface score >= "
                f"{iface_cutoff} and downsample rate {downsample_rate}."
            )

    return selected_points_idx, patch_descs, patch_iface[selected_points_idx]


def multidock(
    source_pt,
    source_pcd,
    source_patch_idxs,
    source_descs,
    target_pt,
    target_pcd,
    target_patch_idxs,
    target_descs,
    binder_align: bool = False,
    ransac_skip: bool = False,
):
    ransac_radius = 1.5
    ransac_iter = 10000
    all_results = []
    all_source_patch = []
    all_source_desc = []
    all_source_idx = []

    target_patch, target_patch_descs, target_patch_idx = get_patch_geo(
        target_pcd,
        target_patch_idxs,
        target_pt,
        target_descs,
        flip_normals=binder_align,
        outward_shift=0.25,
    )

    for pt in source_pt:
        source_patch, source_patch_descs, source_patch_idx = get_patch_geo(
            source_pcd, source_patch_idxs, pt, source_descs, outward_shift=0.25
        )

        if not ransac_skip:
            result = registration_ransac_based_on_feature_matching(
                source=source_patch,
                target=target_patch,
                source_feature=source_patch_descs[0],
                target_feature=target_patch_descs[0],
                max_correspondence_distance=ransac_radius,
                estimation_method=TransformationEstimationPointToPoint(False),
                ransac_n=3,
                checkers=[
                    CorrespondenceCheckerBasedOnEdgeLength(0.9),
                    CorrespondenceCheckerBasedOnDistance(1.0),
                    CorrespondenceCheckerBasedOnNormal(np.pi / 2),
                ],
                criteria=RANSACConvergenceCriteria(ransac_iter, 500),
            )
            init = result.transformation
        else:
            init = np.identity(4)

        result_icp = registration_icp(
            source=source_patch,
            target=target_patch,
            max_correspondence_distance=1.5,
            init=init,
            estimation_method=TransformationEstimationPointToPlane(),
            criteria=ICPConvergenceCriteria(),
        )

        source_patch.transform(result_icp.transformation)
        all_results.append(result_icp)
        all_source_patch.append(source_patch)
        all_source_desc.append(source_patch_descs)
        all_source_idx.append(source_patch_idx)

    return all_results, all_source_patch, all_source_desc, all_source_idx


def transform_patch_coords(pcd, patch_indices, site, T):
    """Transform geodesic patch vertex coordinates for a surface site."""
    pts = np.asarray(pcd.points)[patch_indices[site]]
    return apply_transform(pts, T)

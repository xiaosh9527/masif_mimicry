import os

import numpy as np
import open3d as o3d
from geometry.open3d_import import (
    Feature,
    PointCloud,
    Vector3dVector,
    read_point_cloud,
    read_triangle_mesh,
)
from packaging import version

if version.parse("0.12.0") <= version.parse(o3d.__version__):
    ICPConvergenceCriteria = o3d.pipelines.registration.ICPConvergenceCriteria
elif version.parse("0.6.0") < version.parse(o3d.__version__):
    ICPConvergenceCriteria = o3d.registration.ICPConvergenceCriteria
else:
    ICPConvergenceCriteria = o3d.ICPConvergenceCriteria


def get_features(params: dict, pdb: str, pid: str, source: bool = True, flip_desc: bool = False) -> dict:
    """Load MaSIF features for a PDB id and p1/p2 side."""
    if pid == "p1":
        P = pdb.split("_")[0] + "_" + pdb.split("_")[1]
    else:
        P = pdb.split("_")[0] + "_" + pdb.split("_")[2]

    mask = "seed" if source else "target"

    pdb_fn = os.path.join(params[f"{mask}_pdb_dir"], f"{P}.pdb")
    ply = os.path.join(params[f"{mask}_surf_dir"], f"{P}.ply")
    input_feat_fn = os.path.join(params[f"{mask}_precomp_dir"], pdb, f"{pid}_input_feat.npy")
    indices_fn = os.path.join(params[f"{mask}_precomp_dir"], pdb, f"{pid}_list_indices.npy")
    rho_fn = os.path.join(params[f"{mask}_precomp_dir"], pdb, f"{pid}_rho_wrt_center.npy")
    theta_fn = os.path.join(params[f"{mask}_precomp_dir"], pdb, f"{pid}_theta_wrt_center.npy")
    ilabel_fn = os.path.join(params[f"{mask}_precomp_dir"], pdb, f"{pid}_iface_labels.npy")
    iface_fn = os.path.join(params[f"{mask}_iface_dir"], "pred_" + P + ".npy")
    if flip_desc:
        desc_fn = os.path.join(params[f"{mask}_desc_dir"], pdb, f"{pid}_desc_flipped.npy")
    else:
        desc_fn = os.path.join(params[f"{mask}_desc_dir"], pdb, f"{pid}_desc_straight.npy")

    return dict(
        pdb=pdb_fn,
        mesh=read_triangle_mesh(ply),
        pcd=read_point_cloud(ply),
        rho=np.load(rho_fn),
        theta=np.load(theta_fn),
        desc=np.load(desc_fn),
        input_feat=np.load(input_feat_fn),
        indices=np.load(indices_fn, allow_pickle=True),
        iface=np.load(iface_fn),
        ilabel=np.load(ilabel_fn),
    )


def get_patch_geo(
    pcd,
    patch_coords,
    center,
    descriptors,
    outward_shift=0.25,
    flip_normals=False,
):
    """Return geodesic patch point cloud and descriptors for a surface site."""
    patch_idxs = patch_coords[center]
    patch_pts = np.asarray(pcd.points)[patch_idxs, :]
    patch_nrmls = np.asarray(pcd.normals)[patch_idxs, :]
    patch_pts = patch_pts + outward_shift * patch_nrmls
    if flip_normals:
        patch_nrmls = -patch_nrmls

    patch = PointCloud()
    patch.points = Vector3dVector(patch_pts)
    patch.normals = Vector3dVector(patch_nrmls)
    patch_descs = [Feature(), Feature(), Feature()]
    patch_descs[0].data = descriptors[patch_idxs, :].T

    return patch, patch_descs, patch_idxs

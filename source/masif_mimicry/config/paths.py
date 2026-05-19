import os

from default_config.masif_opts import masif_opts

BENCHMARK_PDB_SUBDIR = "01-benchmark_pdbs"
BENCHMARK_SURF_SUBDIR = "01-benchmark_surfaces"
PRECOMP_12A_SUBDIR = "04b-precomputation_12A/precomputation"
PRECOMP_9A_SUBDIR = "04a-precomputation_9A/precomputation"


def resolve_database_paths(database_dir: str):
    """
    Resolve MaSIF database root and data_preparation directory.

    database_dir may be either the database root (e.g. TED .../output) or
    .../output/data_preparation directly.
    """
    database_dir = os.path.normpath(database_dir)
    if os.path.basename(database_dir) == "data_preparation":
        return os.path.dirname(database_dir), database_dir
    return database_dir, os.path.join(database_dir, "data_preparation")


def set_params(*, database_dir: str, target_preprocess_dir: str, masif_app: str = "ppi_search") -> dict:
    """Set path parameters for mimicry search (seeds + target preprocess trees)."""
    params = {}
    db_root, db_prep = resolve_database_paths(database_dir)
    params["database_root"] = db_root
    params["database_preparation_dir"] = db_prep
    params["top_seed_dir"] = db_prep
    params["masif_target_root"] = os.path.normpath(target_preprocess_dir)
    params["out_dir_template"] = "tmp/{}/"

    params["seed_surf_dir"] = os.path.join(db_prep, BENCHMARK_SURF_SUBDIR)
    params["seed_pdb_dir"] = os.path.join(db_prep, BENCHMARK_PDB_SUBDIR)
    params["seed_precomp_dir"] = os.path.join(db_prep, PRECOMP_12A_SUBDIR)
    params["seed_iface_dir"] = os.path.join(db_root, masif_opts["site"]["out_pred_dir"])
    params["seed_ply_iface_dir"] = os.path.join(db_root, masif_opts["site"]["out_surf_dir"])
    params["seed_desc_dir"] = os.path.join(db_root, masif_opts["ppi_search"]["desc_dir"])

    params["target_surf_dir"] = os.path.join(params["masif_target_root"], masif_opts["ply_chain_dir"])
    params["target_iface_dir"] = os.path.join(params["masif_target_root"], masif_opts["site"]["out_pred_dir"])
    params["target_ply_iface_dir"] = os.path.join(params["masif_target_root"], masif_opts["site"]["out_surf_dir"])
    params["target_pdb_dir"] = os.path.join(params["masif_target_root"], masif_opts["pdb_chain_dir"])
    params["target_desc_dir"] = os.path.join(params["masif_target_root"], masif_opts["ppi_search"]["desc_dir"])
    params["target_precomp_dir"] = os.path.join(
        params["masif_target_root"], masif_opts[masif_app]["masif_precomputation_dir"]
    )

    return params

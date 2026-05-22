"""CLI entry point: python -m masif_mimicry.search"""

import argparse

from masif_mimicry.search.workflow import run_search


def create_parser():
    p = argparse.ArgumentParser("Simplified MaSIF mimicry search")

    p.add_argument(
        "--database_dir",
        type=str,
        required=True,
        help="MaSIF database root (e.g. TED_domainome/output); "
        "seed PDBs live under data_preparation/01-benchmark_pdbs/",
    )
    p.add_argument(
        "--target_preprocess_dir",
        type=str,
        required=True,
        help="MaSIF target preprocess directory (same layout for the target protein)",
    )

    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--seed_pdb", type=str, help="Single seed PDB (format: PDB_ppi_chain or PDB_ppi)")
    group.add_argument("--split_seed_list", type=str, help="File containing one seed PDB id per line")

    p.add_argument(
        "--target_run_dir",
        type=str,
        default=None,
        help="Prepared target run directory (target_sites.json + target_vert/). "
        "When set, target site selection flags are not allowed.",
    )

    p.add_argument(
        "--target_pdb",
        type=str,
        default=None,
        help="Target PDB identifier (legacy mode without --target_run_dir)",
    )
    p.add_argument(
        "--target_ppi_id",
        choices=["p1", "p2"],
        default="p1",
        help="ppi side of the target (legacy mode only; read from manifest if --target_run_dir)",
    )
    p.add_argument(
        "--target_chain",
        type=str,
        default=None,
        help="Target chain id (legacy mode only; read from manifest if --target_run_dir)",
    )

    p.add_argument(
        "--target_residue",
        type=int,
        default=None,
        help="[define_target_sites only] Target residue number",
    )
    p.add_argument(
        "--target_atom",
        type=str,
        default=None,
        help="[define_target_sites only] Target atom name",
    )
    p.add_argument(
        "--num_points",
        type=int,
        default=None,
        help="[define_target_sites only] Number of sites (residue mode) or minimum (exhaustive)",
    )

    p.add_argument(
        "--downsample",
        type=int,
        default=1,
        help="Downsample rate for seed patch selection",
    )
    p.add_argument(
        "--top_iface_percent",
        type=float,
        default=0.0,
        help="Top percentile of seed patches by interface score",
    )
    p.add_argument("--iface_cutoff", type=float, default=0.0, help="Interface score cutoff for seed patches")
    p.add_argument(
        "--interface_only",
        action="store_true",
        help="Search only seed interfaces (requires seed complexes)",
    )
    p.add_argument("--desc_dist_cutoff", type=float, default=1.5, help="Descriptor distance cutoff (filtering)")
    p.add_argument("--desc_dist_score_cutoff", type=float, default=0.25, help="Post-alignment score cutoff")

    p.add_argument("--output_dir", type=str, default=".", help="Output directory (legacy mode only)")
    p.add_argument("--output_postfix", type=str, default="", help="Subfolder postfix (legacy mode only)")
    p.add_argument("--count_clashes", action="store_true", help="Compute clashes (slower)")
    p.add_argument("--ca_clash_threshold", type=float, default=100, help="CA clash threshold")
    p.add_argument("--heavy_atom_clash_threshold", type=float, default=100, help="Heavy-atom clash threshold")
    p.add_argument(
        "--compute_source_residues",
        action="store_true",
        help="Map patch centers to residues via surf2atom (DSSP; slower)",
    )

    return p


def main(argv=None):
    parser = create_parser()
    run_search(parser.parse_args(argv))


if __name__ == "__main__":
    main()

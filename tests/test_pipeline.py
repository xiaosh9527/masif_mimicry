from masif_mimicry.backends.render import render_execution_script
from masif_mimicry.config import JobConfig
from masif_mimicry.pipeline import build_commands


def _job(proteins=None):
    return JobConfig.model_validate(
        {
            "name": "test-run",
            "target": {
                "structure": "/workspace/masif_mimicry/target.pdb",
                "pdb_id": "1abc",
                "query_chain": "A",
                "receptor_chain": "B",
                "focus": {"residue": 10, "atom": "CA"},
                "ligand": "B:LIG",
            },
            "database": {"proteins": proteins or ["seed_A"]},
            "runtime": {"backend": "local", "suite_root": "/workspace"},
            "output_root": "/workspace/masif_mimicry/runs",
        }
    )


def test_commands_use_explicit_python_runtimes():
    commands = "\n".join(build_commands(_job()))
    assert "/usr/local/bin/masif-python" in commands
    assert "/opt/conda/envs/bonsai/bin/python" in commands
    assert "conda activate" not in commands


def test_target_cache_is_shared_between_runs():
    commands = "\n".join(build_commands(_job()))
    assert "/workspace/masif_mimicry/data/target" in commands
    assert "--target_db_root /workspace/masif_mimicry" in commands
    assert "--target_db target" in commands
    assert "/runs/test-run/data/target" in commands
    assert "ln -s ../../../data/target" in commands


def test_single_protein_uses_seed_pdb():
    commands = "\n".join(build_commands(_job()))
    assert "--seed_pdb seed_A" in commands


def test_multiple_proteins_use_seed_list():
    commands = "\n".join(build_commands(_job(["seed_A", "seed_B"])))
    assert "--split_seed_list" in commands


def test_local_execution_script_is_strict_bash():
    script = render_execution_script(_job(), ["echo ok"])
    assert script.startswith("#!/bin/bash")
    assert "set -euo pipefail" in script

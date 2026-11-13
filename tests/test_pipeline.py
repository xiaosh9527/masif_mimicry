import json
from pathlib import Path

import pytest

from masif_mimicry.backends.render import render_execution_script
from masif_mimicry.config import JobConfig
from masif_mimicry.pipeline import build_commands, database_proteins, prepare_run


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


def _array_job(tmp_path: Path, proteins: list[str], tasks: int = 2) -> JobConfig:
    suite = tmp_path / "masif_seed"
    target = suite / "masif_mimicry" / "inputs" / "target.pdb"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("", encoding="utf-8")
    return JobConfig.model_validate(
        {
            "name": "proteome-scan",
            "target": {
                "structure": str(target),
                "pdb_id": "1abc",
                "query_chain": "A",
            },
            "database": {"proteins": proteins},
            "postprocessing": {"enabled": False},
            "runtime": {
                "backend": "slurm",
                "suite_root": str(suite),
                "image": str(suite / "masif_mimicry.sif"),
                "slurm": {
                    "nodes": 1,
                    "qos": "serial",
                    "ntasks_per_node": 1,
                    "job_name": "benchmark",
                    "array_tasks": tasks,
                    "cpus": 2,
                    "memory": "10g",
                    "time": "24:00:00",
                    "output": "logs/%x-%A_%a.out",
                },
            },
            "output_root": str(suite / "masif_mimicry" / "runs"),
        }
    )


def test_database_proteins_can_come_from_manifest(tmp_path: Path):
    manifest = tmp_path / "proteins.txt"
    manifest.write_text("# stable seed list\nseed_A\n\nseed_B\n", encoding="utf-8")
    job = _job()
    data = job.model_dump()
    data["database"] = {"protein_list": manifest}

    assert database_proteins(JobConfig.model_validate(data)) == ["seed_A", "seed_B"]


def test_slurm_array_header_matches_requested_resources(tmp_path: Path):
    script = render_execution_script(_array_job(tmp_path, ["a", "b"]), ["echo ok"])

    for directive in (
        "#!/bin/bash -l",
        "#SBATCH --nodes=1",
        "#SBATCH --qos=serial",
        "#SBATCH --ntasks-per-node=1",
        "#SBATCH --job-name=benchmark",
        "#SBATCH --array=0-1",
        "#SBATCH --cpus-per-task=2",
        "#SBATCH --mem=10g",
        "#SBATCH --time=24:00:00",
        "#SBATCH --output=logs/%x-%A_%a.out",
    ):
        assert directive in script


def test_apptainer_additional_bind_is_rendered(tmp_path: Path):
    job = _array_job(tmp_path, ["a", "b"])
    external = tmp_path / "external-database"
    external.mkdir()
    job.runtime.binds = [external]

    script = render_execution_script(job, ["echo ok"])

    assert f"--bind {external}:{external}" in script


def test_prepare_run_writes_deterministic_restartable_shards(tmp_path: Path):
    proteins = ["seed_A", "seed_B", "seed_C", "seed_D", "seed_E"]
    job = _array_job(tmp_path, proteins)

    plan = prepare_run(job)

    shard_dir = plan.run_dir / "inputs" / "seed-shards"
    assert (shard_dir / "0.txt").read_text(encoding="utf-8") == "seed_A\nseed_C\nseed_E\n"
    assert (shard_dir / "1.txt").read_text(encoding="utf-8") == "seed_B\nseed_D\n"
    array_plan = json.loads(
        (plan.run_dir / "inputs" / "array-plan.json").read_text(encoding="utf-8")
    )
    assert array_plan["protein_count"] == 5
    assert array_plan["array_tasks"] == 2
    script = plan.script_path.read_text(encoding="utf-8")
    assert 'SLURM_ARRAY_TASK_ID' in script
    assert 'Shard $MASIF_ARRAY_TASK_ID already complete; skipping.' in script
    assert 'touch "$MASIF_SHARD_DONE"' in script
    assert 'results/shards' not in script
    assert '--output_dir /workspace/masif_mimicry/runs/proteome-scan/results' in script
    assert '--flat_output' in script
    assert 'tmp/shards' in script


def test_completed_array_cannot_be_reassigned(tmp_path: Path):
    first = _array_job(tmp_path, ["seed_A", "seed_B"])
    plan = prepare_run(first)
    marker = plan.run_dir / "state" / "shards" / "0.complete"
    marker.write_text("", encoding="utf-8")
    changed = _array_job(tmp_path, ["seed_A", "seed_C"])

    with pytest.raises(ValueError, match="different protein list or array size"):
        prepare_run(changed)

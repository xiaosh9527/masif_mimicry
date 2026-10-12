from __future__ import annotations

import hashlib
import json
import os
import platform
import shlex
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .backends.render import command_path, render_execution_script
from .config import JobConfig, dump_job


def _q(value: object) -> str:
    return shlex.quote(str(value))


@dataclass(frozen=True)
class RunPlan:
    run_dir: Path
    script_path: Path
    commands: tuple[str, ...]


def _project_root(job: JobConfig) -> Path:
    return job.runtime.suite_root / "masif_mimicry"


def _python_runtimes() -> tuple[str, str]:
    """Return the isolated MaSIF 3.6 and bonsai 3.11 interpreters."""
    return (
        os.getenv("MASIF_PYTHON", "/usr/local/bin/masif-python"),
        os.getenv("BONSAI_PYTHON", "/opt/conda/envs/bonsai/bin/python"),
    )


def database_proteins(job: JobConfig) -> list[str]:
    """Return the deterministic seed list from inline config or a manifest file."""
    if job.database.proteins:
        proteins = list(job.database.proteins)
    else:
        assert job.database.protein_list is not None
        try:
            lines = job.database.protein_list.read_text(encoding="utf-8").splitlines()
        except OSError as exc:
            raise ValueError(
                f"cannot read database.protein_list: {job.database.protein_list}"
            ) from exc
        proteins = [
            line.strip()
            for line in lines
            if line.strip() and not line.lstrip().startswith("#")
        ]
    if not proteins:
        raise ValueError("database protein source contains no identifiers")
    if len(set(proteins)) != len(proteins):
        raise ValueError("database protein source contains duplicate identifiers")
    return proteins


def _replace_shell_tokens(command: str) -> str:
    return command.replace("__MASIF_SEED_SHARD__", '"$MASIF_SEED_SHARD"')


def build_commands(job: JobConfig) -> list[str]:
    proteins = database_proteins(job)
    array_tasks = job.runtime.slurm.array_tasks if job.runtime.backend == "slurm" else None
    array_scan = array_tasks is not None
    project_root = _project_root(job)
    run_dir = job.run_dir
    target_dir = project_root / "data" / "target"
    target_link = run_dir / "data" / "target"
    result_dir = run_dir / "results"
    tmp_dir = run_dir / "tmp"
    database_root = job.database.root or (job.runtime.suite_root / "masif")

    c_project = command_path(job, project_root)
    c_run = command_path(job, run_dir)
    c_target_structure = command_path(job, job.target.structure)
    c_target_dir = command_path(job, target_dir)
    c_target_link = command_path(job, target_link)
    c_target_link_value = os.path.relpath(c_target_dir, start=Path(c_target_link).parent)
    c_result_dir = command_path(job, result_dir)
    c_tmp = command_path(job, tmp_dir)
    c_database_root = command_path(job, database_root)
    masif_python, bonsai_python = _python_runtimes()

    source_paths = [
        Path(c_database_root) / "source",
        Path("/workspace/masif_seed_search/source")
        if job.runtime.backend != "local"
        else job.runtime.suite_root / "masif_seed_search" / "source",
        Path(c_project) / "source",
    ]
    pythonpath = ":".join(str(path) for path in source_paths)
    link_command = (
        f"if [[ ! -e {_q(c_target_link)} ]]; then "
        f"ln -s {_q(c_target_link_value)} {_q(c_target_link)} 2>/dev/null || "
        f"[[ -L {_q(c_target_link)} ]]; fi"
    )
    commands = [
        (
            f"test -x {_q(masif_python)} || "
            f"{{ echo 'MaSIF Python not found: {masif_python}' >&2; exit 1; }}"
        ),
        (
            f"test -x {_q(bonsai_python)} || "
            f"{{ echo 'Bonsai Python not found: {bonsai_python}' >&2; exit 1; }}"
        ),
        f"mkdir -p {_q(c_target_dir)} {_q(c_result_dir)} {_q(c_tmp)} {_q(Path(c_target_link).parent)}",
        link_command,
        f"export MASIF_PYTHON={_q(masif_python)}",
        f"export BONSAI_PYTHON={_q(bonsai_python)}",
        f"export PYTHONPATH={_q(pythonpath)}:${{PYTHONPATH:-}}",
    ]
    if array_scan:
        c_shard_dir = command_path(job, run_dir / "inputs" / "seed-shards")
        c_state_dir = command_path(job, run_dir / "state" / "shards")
        c_shard_tmp = command_path(job, tmp_dir / "shards")
        commands.extend(
            [
                'export MASIF_ARRAY_TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"',
                f"printf -v MASIF_SEED_SHARD '%s/%s.txt' {_q(c_shard_dir)} \"$MASIF_ARRAY_TASK_ID\"",
                f"printf -v MASIF_SHARD_TMP_DIR '%s/%s' {_q(c_shard_tmp)} \"$MASIF_ARRAY_TASK_ID\"",
                f"printf -v MASIF_SHARD_DONE '%s/%s.complete' {_q(c_state_dir)} \"$MASIF_ARRAY_TASK_ID\"",
                'if [[ -f "$MASIF_SHARD_DONE" ]]; then echo "Shard $MASIF_ARRAY_TASK_ID already complete; skipping."; exit 0; fi',
                'test -s "$MASIF_SEED_SHARD" || { echo "Missing or empty shard: $MASIF_SEED_SHARD" >&2; exit 1; }',
                'mkdir -p "$MASIF_SHARD_TMP_DIR"',
                'export TMPDIR="$MASIF_SHARD_TMP_DIR"',
            ]
        )
    else:
        commands.append(f"export TMPDIR={_q(c_tmp)}")
    commands.append(
        f"bash {_q(Path(c_project) / 'preprocess_pdb.sh')} "
        f"{_q(c_target_structure)} {_q(job.target.processed_id)} -o {_q(c_target_dir)}"
    )

    params = job.search.resolved()
    search = [
        masif_python,
        "-u",
        str(Path(c_project) / "source" / "masif_mimicry_search.py"),
        "--seed_db_root",
        c_database_root,
        "--seed_db",
        job.database.name,
        "--target_db_root",
        c_project,
        "--target_db",
        "target",
        "--target_pdb",
        job.target.processed_id,
        "--target_chain",
        job.target.query_chain,
        "--downsample",
        params["downsample"],
        "--desc_dist_cutoff",
        params["descriptor_distance_cutoff"],
        "--desc_dist_score_cutoff",
        params["mimicry_score_cutoff"],
        "--iface_cutoff",
        params["interface_cutoff"],
        "--top_iface_percent",
        params["top_interface_percent"],
        "--num_points",
        params["number_of_points"],
        "--ca_clash_threshold",
        params["ca_clash_threshold"],
        "--heavy_atom_clash_threshold",
        params["heavy_atom_clash_threshold"],
        "--output_dir",
        c_result_dir,
        "--output_postfix",
        job.name,
    ]
    if array_scan:
        search.extend(["--split_seed_list", "__MASIF_SEED_SHARD__", "--flat_output"])
    elif len(proteins) == 1:
        search.extend(["--seed_pdb", proteins[0]])
    else:
        seed_list = run_dir / "inputs" / "seed-proteins.txt"
        search.extend(["--split_seed_list", command_path(job, seed_list)])
    if job.target.focus:
        search.extend(
            [
                "--target_residue",
                job.target.focus.residue,
                "--target_atom",
                job.target.focus.atom,
            ]
        )
    if params["count_clashes"]:
        search.append("--count_clashes")
    if params["interface_only"]:
        search.append("--interface_only")
    commands.append(_replace_shell_tokens(" ".join(_q(item) for item in search)))

    if job.postprocessing.enabled:
        ligand = job.target.ligand.replace(":", "_") if job.target.ligand else None
        if ligand is None:
            commands.append(
                "echo 'Postprocessing skipped: target.ligand is not configured.' >&2"
            )
        else:
            search_root = Path(c_result_dir) / f"{job.target.processed_id}_{job.name}"
            postprocess_dir = Path(c_run) / "postprocess"
            trunc_dir = Path(c_run) / f"truncate-{job.postprocessing.truncation_length}"
            postprocess_script = (
                Path(c_project) / "scripts" / "bonsai_scripts" / "python" / "postprocess_masif_mimicry.py"
            )
            truncate_script = (
                Path(c_project) / "scripts" / "bonsai_scripts" / "python" / "EvoEF2_truncate.py"
            )
            rescore_script = (
                Path(c_project) / "scripts" / "bonsai_scripts" / "python" / "proc_trunc_masif_mimicry.py"
            )
            postprocessed_csv = postprocess_dir / "postprocessed_scores.csv"
            truncated_csv = postprocess_dir / "postprocessed_scores_truncated_results.csv"
            commands.extend(
                [
                    f"mkdir -p {_q(postprocess_dir)} {_q(trunc_dir / 'structures')} {_q(trunc_dir / 'tmp')}",
                    (
                        "pose_paths=(); while IFS= read -r -d '' pose; do pose_paths+=(\"$pose\"); done "
                        f"< <(find {_q(search_root)} -mindepth 2 -type f -name '*.pdb' -print0)"
                    ),
                    (
                        "if (( ${#pose_paths[@]} == 0 )); then "
                        "echo 'No accepted poses were produced; skipping postprocessing.' >&2; exit 0; fi"
                    ),
                    (
                        f"{_q(bonsai_python)} {_q(postprocess_script)} --pdb_paths "
                        f'"${{pose_paths[@]}}" --target_pdb_path {_q(c_target_structure)} '
                        f"--ligand {_q(ligand)} -o {_q(postprocessed_csv)}"
                    ),
                    (
                        f"{_q(bonsai_python)} {_q(truncate_script)} --clean "
                        f"--input_csv {_q(postprocessed_csv)} "
                        f"--length {job.postprocessing.truncation_length} "
                        f"--out_dir {_q(trunc_dir / 'structures')} "
                        f"--temp_dir {_q(trunc_dir / 'tmp')}"
                    ),
                    (
                        f"{_q(bonsai_python)} {_q(rescore_script)} "
                        f"--input {_q(truncated_csv)} --ligand {_q(ligand)} "
                        f"-o {_q(trunc_dir / 'results.csv')}"
                    ),
                ]
            )
    if array_scan:
        commands.append('touch "$MASIF_SHARD_DONE"')
    return commands


def prepare_run(job: JobConfig) -> RunPlan:
    proteins = database_proteins(job)
    array_tasks = job.runtime.slurm.array_tasks if job.runtime.backend == "slurm" else None
    run_dir = job.run_dir
    for relative in ("inputs", "logs", "results", "tmp", "data", "state"):
        (run_dir / relative).mkdir(parents=True, exist_ok=True)
    if array_tasks is not None:
        if array_tasks > len(proteins):
            raise ValueError(
                f"slurm.array_tasks ({array_tasks}) exceeds the number of proteins "
                f"({len(proteins)})"
            )
        shard_dir = run_dir / "inputs" / "seed-shards"
        state_dir = run_dir / "state" / "shards"
        shard_dir.mkdir(parents=True, exist_ok=True)
        state_dir.mkdir(parents=True, exist_ok=True)
        protein_digest = hashlib.sha256(
            ("\n".join(proteins) + "\n").encode("utf-8")
        ).hexdigest()
        array_plan = {
            "format": 1,
            "assignment": "round-robin",
            "protein_count": len(proteins),
            "protein_sha256": protein_digest,
            "array_tasks": array_tasks,
        }
        array_plan_path = run_dir / "inputs" / "array-plan.json"
        complete_markers = list(state_dir.glob("*.complete"))
        if complete_markers:
            try:
                previous_plan = json.loads(array_plan_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                raise ValueError(
                    f"run {job.name!r} has completed shards but no readable array plan; "
                    "use a new run name"
                ) from exc
            if previous_plan != array_plan:
                raise ValueError(
                    f"run {job.name!r} has completed shards for a different protein list "
                    "or array size; use a new run name"
                )
        array_plan_path.write_text(
            json.dumps(array_plan, indent=2) + "\n", encoding="utf-8"
        )
        for shard_id in range(array_tasks):
            shard = proteins[shard_id::array_tasks]
            (shard_dir / f"{shard_id}.txt").write_text(
                "\n".join(shard) + "\n", encoding="utf-8"
            )
    elif len(proteins) > 1:
        (run_dir / "inputs" / "seed-proteins.txt").write_text(
            "\n".join(proteins) + "\n", encoding="utf-8"
        )

    resolved_path = run_dir / "run.yaml"
    dump_job(job, resolved_path)
    commands = build_commands(job)
    script_name = "submit.slurm" if job.runtime.backend == "slurm" else "run.sh"
    script_path = run_dir / script_name
    script_path.write_text(render_execution_script(job, commands), encoding="utf-8")
    script_path.chmod(0o755)
    provenance = {
        "masif_mimicry_interface_version": __version__,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "host": platform.node(),
        "platform": platform.platform(),
        "backend": job.runtime.backend,
        "pid": os.getpid(),
        "database_proteins": len(proteins),
        "slurm_array_tasks": array_tasks,
    }
    (run_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    return RunPlan(run_dir=run_dir, script_path=script_path, commands=tuple(commands))


def execute(plan: RunPlan, job: JobConfig) -> subprocess.CompletedProcess[str]:
    if job.runtime.backend == "slurm":
        command = (
            ["sbatch", plan.script_path.name]
            if job.runtime.submit
            else ["bash", "-n", plan.script_path.name]
        )
        return subprocess.run(command, cwd=plan.run_dir, text=True, check=True)
    else:
        command = [str(plan.script_path)]
    return subprocess.run(command, text=True, check=True)

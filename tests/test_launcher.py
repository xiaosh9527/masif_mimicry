"""Host launcher tests: containers and the scheduler are always inert stubs."""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml


@pytest.fixture
def launcher(tmp_path):
    project = tmp_path / "suite" / "masif_mimicry"
    (project / "bin").mkdir(parents=True)
    script = project / "bin" / "masif-mimicry"
    shutil.copyfile(Path(__file__).parents[1] / "bin" / "masif-mimicry", script)
    image = project / "image.sif"
    image.touch()
    stubs = tmp_path / "commands"
    stubs.mkdir()
    # Isolate PATH so these tests can never find a real sbatch or container runtime.
    for command in ("bash", "realpath", "dirname"):
        (stubs / command).symlink_to(shutil.which(command))
    apptainer = stubs / "apptainer"
    apptainer.write_text(
        f"#!{sys.executable}\n"
        "import json, os, pathlib, subprocess, sys\n"
        "args = sys.argv[1:]\n"
        "if '-c' in args:\n"
        "    i = args.index('-c')\n"
        "    sys.exit(subprocess.call([sys.executable, *args[i:]]))\n"
        "pathlib.Path(os.environ['PREPARE_LOG']).write_text(json.dumps(args))\n"
        "if os.environ.get('PREPARE_FAIL'):\n"
        "    sys.exit(7)\n"
        "run = pathlib.Path(os.environ['RUN_DIR'])\n"
        "run.mkdir(parents=True, exist_ok=True)\n"
        "if not os.environ.get('MISSING_SCRIPT'):\n"
        "    (run / 'submit.slurm').write_text(os.environ['SCRIPT_CONTENT'])\n"
    )
    apptainer.chmod(0o755)
    sbatch = stubs / "sbatch"
    sbatch.write_text(
        f"#!{sys.executable}\n"
        "import json, os, pathlib, sys\n"
        "pathlib.Path(os.environ['SBATCH_LOG']).write_text(\n"
        "    json.dumps({'args': sys.argv[1:], 'cwd': os.getcwd()}))\n"
    )
    sbatch.chmod(0o755)
    run_dir = project / "runs" / "example"
    env = {
        **os.environ,
        "PATH": str(stubs),
        "MIMICRY_IMAGE": str(image),
        "SBATCH_LOG": str(tmp_path / "sbatch.json"),
        "PREPARE_LOG": str(tmp_path / "prepare.json"),
        "RUN_DIR": str(run_dir),
        "SCRIPT_CONTENT": "#!/bin/bash\ntrue\n",
    }
    for key in ("PREPARE_FAIL", "MISSING_SCRIPT"):
        env.pop(key, None)

    def run(submit=None, dry_run=False, *, scheduler=True, **overrides):
        runtime = {
            "suite_root": str(project.parent),
            "image": str(image),
            "backend": "slurm",
            "binds": [str(tmp_path / "extra bind")],
        }
        if submit is not None:
            runtime["submit"] = submit
        job = project / "job.yaml"
        job.write_text(yaml.safe_dump({"name": "example", "runtime": runtime}))
        if not scheduler:
            sbatch.unlink()
        args = [str(stubs / "bash"), str(script), "run"]
        if dry_run:
            args.append("--dry-run")
        return subprocess.run(
            [*args, str(job)],
            env={**env, **overrides},
            text=True,
            capture_output=True,
            check=False,
        )

    return run, env


@pytest.mark.parametrize("submit", [None, False, True, "false", "true"])
@pytest.mark.parametrize("dry_run", [False, True])
def test_host_slurm_requires_opt_in(launcher, submit, dry_run):
    run, env = launcher
    result = run(submit, dry_run)
    assert result.returncode == 0, result.stderr
    log = Path(env["SBATCH_LOG"])
    should_submit = submit in (True, "true") and not dry_run
    assert log.exists() == should_submit
    if should_submit:
        assert json.loads(log.read_text()) == {"args": ["submit.slurm"], "cwd": env["RUN_DIR"]}
    else:
        assert "Prepared without submission:" in result.stdout
    preparation = json.loads(Path(env["PREPARE_LOG"]).read_text())
    assert "--dry-run" in preparation
    bind = str(Path(env["PREPARE_LOG"]).parent / "extra bind")
    assert f"{bind}:{bind}" in preparation


@pytest.mark.parametrize("submit,dry_run", [(None, False), (False, False), (True, True)])
def test_preparation_does_not_require_sbatch(launcher, submit, dry_run):
    run, env = launcher
    result = run(submit, dry_run, scheduler=False)
    assert result.returncode == 0, result.stderr
    assert not Path(env["SBATCH_LOG"]).exists()


@pytest.mark.parametrize("submit", [False, True])
@pytest.mark.parametrize("failure", ["syntax", "missing", "preparation"])
def test_invalid_preparation_never_submits(launcher, submit, failure):
    run, env = launcher
    overrides = {
        "syntax": {"SCRIPT_CONTENT": "if then\n"},
        "missing": {"MISSING_SCRIPT": "1"},
        "preparation": {"PREPARE_FAIL": "1"},
    }[failure]
    result = run(submit, **overrides)
    assert result.returncode != 0
    assert not Path(env["SBATCH_LOG"]).exists()


def test_submission_requires_host_sbatch(launcher):
    run, env = launcher
    result = run(True, scheduler=False)
    assert result.returncode != 0
    assert "Slurm sbatch was not found on the host." in result.stderr
    assert not Path(env["SBATCH_LOG"]).exists()

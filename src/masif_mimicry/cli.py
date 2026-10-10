from __future__ import annotations

import json
import shutil
import subprocess
import sys
import webbrowser
from pathlib import Path
from typing import Annotated, Optional

import typer
from pydantic import ValidationError

from . import __version__
from .checkpoint import inspect_checkpoint
from .config import load_job
from .doctor import all_ok, run_checks
from .pipeline import build_commands, database_proteins, execute, prepare_run
from .results import build_report

app = typer.Typer(
    no_args_is_help=True,
    invoke_without_command=True,
    rich_markup_mode="markdown",
    help="Run MaSIF surface-mimicry workflows without editing shell scripts.",
)


def _load(path: Path):
    try:
        return load_job(path)
    except (OSError, ValueError, ValidationError) as exc:
        typer.secho(f"Configuration error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(2) from exc


@app.callback()
def main(version: Annotated[bool, typer.Option("--version", help="Show the version.")] = False):
    if version:
        typer.echo(__version__)
        raise typer.Exit()


@app.command("init")
def init_project(
    directory: Annotated[Path, typer.Argument(help="Directory for the new experiment.")],
    force: Annotated[bool, typer.Option(help="Replace an existing example job file.")] = False,
):
    directory = directory.expanduser().resolve()
    target = directory / "job.yaml"
    if target.exists() and not force:
        typer.secho(f"Refusing to replace {target}; use --force.", fg=typer.colors.RED, err=True)
        raise typer.Exit(2)
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "inputs").mkdir(exist_ok=True)
    example = Path(__file__).parent / "templates" / "job.yaml"
    shutil.copyfile(example, target)
    typer.echo(f"Created {target}")
    typer.echo("Put the target PDB in inputs/, edit job.yaml, then run: masif-mimicry validate job.yaml")


@app.command()
def validate(job_file: Annotated[Path, typer.Argument(exists=True, dir_okay=False)]):
    job = _load(job_file)
    try:
        proteins = database_proteins(job)
    except ValueError as exc:
        typer.secho(f"Configuration error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(2) from exc
    typer.secho(f"Valid job: {job.name}", fg=typer.colors.GREEN)
    typer.echo(f"Target: {job.target.processed_id}")
    typer.echo(f"Database proteins: {len(proteins)}")
    typer.echo(f"Backend: {job.runtime.backend}")


@app.command()
def doctor(job_file: Annotated[Optional[Path], typer.Argument()] = None):
    job = _load(job_file) if job_file else None
    checks = run_checks(job)
    for check in checks:
        typer.secho("✓" if check.ok else "✗", fg=typer.colors.GREEN if check.ok else typer.colors.RED, nl=False)
        typer.echo(f" {check.name}: {check.detail}")
    if not all_ok(checks):
        raise typer.Exit(1)


@app.command()
def checkpoint(path: Annotated[Path, typer.Argument(help="TensorFlow checkpoint prefix or file.")]):
    """Reject model checkpoints containing NaN or infinite weights."""
    try:
        report = inspect_checkpoint(path)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    typer.echo(f"Checkpoint: {report.prefix}")
    typer.echo(f"Tensors: {report.tensors}")
    typer.echo(f"Floating-point values: {report.floating_values}")
    typer.echo(f"Non-finite values: {report.nonfinite_values}")
    if not report.ok:
        typer.secho("Checkpoint is not usable.", fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    typer.secho("Checkpoint is finite.", fg=typer.colors.GREEN)


@app.command()
def plan(job_file: Annotated[Path, typer.Argument(exists=True, dir_okay=False)]):
    job = _load(job_file)
    typer.echo(f"Run directory: {job.run_dir}")
    try:
        commands = build_commands(job)
    except ValueError as exc:
        typer.secho(f"Cannot build workflow: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(2) from exc
    for number, command in enumerate(commands, start=1):
        typer.echo(f"{number}. {command}")


@app.command()
def run(
    job_file: Annotated[Path, typer.Argument(exists=True, dir_okay=False)],
    dry_run: Annotated[bool, typer.Option("--dry-run", help="Generate files without executing.")] = False,
):
    job = _load(job_file)
    checks = run_checks(job)
    failed = [check for check in checks if not check.ok]
    if failed and not dry_run:
        for check in failed:
            typer.secho(f"✗ {check.name}: {check.detail}", fg=typer.colors.RED, err=True)
        typer.echo("Fix these checks or use --dry-run to inspect generated files.", err=True)
        raise typer.Exit(1)
    try:
        plan_data = prepare_run(job)
    except ValueError as exc:
        typer.secho(f"Cannot prepare workflow: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(2) from exc
    typer.echo(f"Prepared {plan_data.run_dir}")
    typer.echo(f"Execution script: {plan_data.script_path}")
    if dry_run:
        return
    try:
        result = execute(plan_data, job)
    except subprocess.CalledProcessError as exc:
        typer.secho(
            f"Workflow command failed with exit code {exc.returncode}. See {plan_data.run_dir / 'logs'}.",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(exc.returncode) from exc
    if result.stdout and result.stdout.strip():
        typer.echo(result.stdout.strip())
    if job.runtime.backend == "slurm" and not job.runtime.submit:
        typer.echo("Slurm script validated but not submitted. Set runtime.submit: true to submit it.")


@app.command()
def results(
    run_dir: Annotated[Path, typer.Argument(exists=True, file_okay=False)],
    open_browser: Annotated[bool, typer.Option("--open", help="Open the generated report.")] = False,
):
    report = build_report(run_dir.expanduser().resolve())
    typer.echo(f"Wrote {report}")
    if open_browser:
        webbrowser.open(report.as_uri())


@app.command()
def status(run_dir: Annotated[Path, typer.Argument(exists=True, file_okay=False)]):
    """Show the files and current stage of a prepared run."""
    run_dir = run_dir.expanduser().resolve()
    stages = [
        ("Prepared", run_dir / "run.yaml"),
        ("Target preprocessed", run_dir / "data" / "target" / "targets"),
        ("Search results", run_dir / "results"),
        ("Postprocessed", run_dir / "postprocess" / "postprocessed_scores.csv"),
    ]
    for label, path in stages:
        complete = path.exists() and (not path.is_dir() or any(path.iterdir()))
        typer.secho("✓" if complete else "○", fg=typer.colors.GREEN if complete else None, nl=False)
        typer.echo(f" {label}")
    shard_dir = run_dir / "inputs" / "seed-shards"
    if shard_dir.is_dir():
        total_shards = len(list(shard_dir.glob("*.txt")))
        complete_shards = len(list((run_dir / "state" / "shards").glob("*.complete")))
        typer.echo(f"Array shards: {complete_shards}/{total_shards} complete")
    provenance_path = run_dir / "provenance.json"
    if provenance_path.exists():
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        typer.echo(f"Backend: {provenance.get('backend', 'unknown')}")


@app.command()
def ui(
    no_open: Annotated[
        bool, typer.Option("--no-open", help="Do not open a browser automatically.")
    ] = False,
):
    """Open the optional browser interface."""
    try:
        import streamlit  # noqa: F401
    except ImportError as exc:
        typer.secho(
            "The browser interface is optional. Install it with: pip install -e '.[ui]'",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(1) from exc
    module_path = Path(__file__).with_name("ui.py")
    command = [sys.executable, "-m", "streamlit", "run", str(module_path)]
    if no_open:
        command.extend(["--server.headless", "true"])
    raise typer.Exit(subprocess.call(command))


if __name__ == "__main__":
    app()

from __future__ import annotations

import shlex
from pathlib import Path

from masif_mimicry.config import JobConfig


def _slurm_header(job: JobConfig) -> str:
    resources = job.runtime.slurm
    lines = [
        "#!/bin/bash -l",
        f"#SBATCH --job-name={job.name}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={resources.cpus}",
        f"#SBATCH --mem={resources.memory}",
        f"#SBATCH --time={resources.time}",
        f"#SBATCH --output={job.run_dir}/logs/slurm-%j.out",
    ]
    if resources.partition:
        lines.append(f"#SBATCH --partition={resources.partition}")
    return "\n".join(lines)


def _container_path(path: Path, suite_root: Path) -> Path:
    try:
        return Path("/workspace") / path.resolve().relative_to(suite_root.resolve())
    except ValueError as exc:
        raise ValueError(
            f"{path} is outside runtime.suite_root ({suite_root}); "
            "container backends can only access files below suite_root"
        ) from exc


def render_execution_script(job: JobConfig, commands: list[str]) -> str:
    body = "set -euo pipefail\n" + "\n".join(commands) + "\n"
    if job.runtime.backend == "local":
        return "#!/bin/bash\n" + body

    assert job.runtime.image is not None
    suite = job.runtime.suite_root.resolve()
    image = job.runtime.image.resolve()
    inner = body
    wrapper = (
        f"apptainer exec --bind {shlex.quote(str(suite))}:/workspace "
        f"--pwd /workspace {shlex.quote(str(image))} /bin/bash -lc {shlex.quote(inner)}"
    )
    if job.runtime.backend == "apptainer":
        return "#!/bin/bash\nset -euo pipefail\n" + wrapper + "\n"
    return _slurm_header(job) + "\n\nset -euo pipefail\n" + wrapper + "\n"


def command_path(job: JobConfig, path: Path) -> str:
    if job.runtime.backend == "local":
        return str(path.resolve())
    return str(_container_path(path, job.runtime.suite_root))


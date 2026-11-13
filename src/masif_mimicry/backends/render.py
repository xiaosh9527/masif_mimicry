from __future__ import annotations

import shlex
from pathlib import Path

from masif_mimicry.config import JobConfig


def _slurm_header(job: JobConfig) -> str:
    resources = job.runtime.slurm
    job_name = resources.job_name or job.name
    output = resources.output or (
        "logs/%x-%A_%a.out" if resources.array_tasks is not None else "logs/%x-%j.out"
    )
    lines = [
        "#!/bin/bash -l",
        f"#SBATCH --nodes={resources.nodes}",
        f"#SBATCH --ntasks-per-node={resources.ntasks_per_node}",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --cpus-per-task={resources.cpus}",
        f"#SBATCH --mem={resources.memory}",
        f"#SBATCH --time={resources.time}",
        f"#SBATCH --output={output}",
    ]
    if resources.qos:
        lines.insert(2, f"#SBATCH --qos={resources.qos}")
    if resources.partition:
        lines.append(f"#SBATCH --partition={resources.partition}")
    if resources.array_tasks is not None:
        array = f"0-{resources.array_tasks - 1}"
        if resources.array_max_parallel is not None:
            array += f"%{resources.array_max_parallel}"
        lines.append(f"#SBATCH --array={array}")
    return "\n".join(lines)


def _container_path(path: Path, suite_root: Path, binds: list[Path]) -> Path:
    try:
        return Path("/workspace") / path.absolute().relative_to(suite_root.absolute())
    except ValueError:
        resolved = path.resolve()
        if any(resolved.is_relative_to(bind.resolve()) for bind in binds):
            return resolved
        raise ValueError(
            f"{path} is outside runtime.suite_root ({suite_root}) and runtime.binds"
        ) from None


def render_execution_script(job: JobConfig, commands: list[str]) -> str:
    body = "set -euo pipefail\n" + "\n".join(commands) + "\n"
    if job.runtime.backend == "local":
        return "#!/bin/bash\n" + body

    assert job.runtime.image is not None
    suite = job.runtime.suite_root.resolve()
    image = job.runtime.image.resolve()
    inner = body
    bind_specs = [f"{suite}:/workspace"]
    bind_specs.extend(f"{path.resolve()}:{path.resolve()}" for path in job.runtime.binds)
    bind_arguments = " ".join(
        f"--bind {shlex.quote(bind_spec)}" for bind_spec in bind_specs
    )
    wrapper = (
        f"apptainer exec {bind_arguments} "
        f"--pwd /workspace {shlex.quote(str(image))} /bin/bash -lc {shlex.quote(inner)}"
    )
    if job.runtime.backend == "apptainer":
        return "#!/bin/bash\nset -euo pipefail\n" + wrapper + "\n"
    return _slurm_header(job) + "\n\nset -euo pipefail\n" + wrapper + "\n"


def command_path(job: JobConfig, path: Path) -> str:
    if job.runtime.backend == "local":
        return str(path.resolve())
    return str(_container_path(path, job.runtime.suite_root, job.runtime.binds))

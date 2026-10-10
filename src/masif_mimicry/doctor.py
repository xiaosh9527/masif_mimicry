from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from .config import JobConfig


@dataclass(frozen=True)
class Check:
    name: str
    ok: bool
    detail: str


def _path_check(name: str, path: Path, kind: str = "path") -> Check:
    return Check(name, path.exists(), f"{path} ({kind})")


def _python_check(name: str, path: str, expected: tuple[int, int]) -> Check:
    executable = Path(path)
    if not executable.is_file():
        return Check(name, False, f"{path} (not found)")
    probe = "import sys; print('.'.join(map(str, sys.version_info[:2])))"
    try:
        result = subprocess.run(
            [path, "-c", probe], check=True, capture_output=True, text=True
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        return Check(name, False, f"{path} ({exc})")
    actual = result.stdout.strip()
    wanted = ".".join(map(str, expected))
    return Check(name, actual == wanted, f"{path} (Python {actual}; expected {wanted})")


def run_checks(job: JobConfig | None = None) -> list[Check]:
    checks = [
        Check("CLI Python", sys.version_info[:2] >= (3, 9), sys.executable),
        Check("Git", shutil.which("git") is not None, shutil.which("git") or "not found"),
    ]
    if job is None:
        return checks
    checks.append(_path_check("Target structure", job.target.structure, "input PDB"))
    if job.database.protein_list is not None:
        checks.append(
            _path_check(
                "Database protein list", job.database.protein_list, "seed manifest"
            )
        )
    checks.append(_path_check("Suite root", job.runtime.suite_root, "suite directory"))
    checks.append(
        _path_check(
            "MaSIF source", job.runtime.suite_root / "masif" / "source", "dependency"
        )
    )
    if job.runtime.backend == "local":
        checks.extend(
            [
                _python_check(
                    "MaSIF runtime",
                    os.getenv("MASIF_PYTHON", "/usr/local/bin/masif-python"),
                    (3, 6),
                ),
                _python_check(
                    "Bonsai runtime",
                    os.getenv("BONSAI_PYTHON", "/opt/conda/envs/bonsai/bin/python"),
                    (3, 11),
                ),
            ]
        )
    checks.append(
        _path_check(
            "MaSIF seed source",
            job.runtime.suite_root / "masif_seed_search" / "source",
            "dependency",
        )
    )
    if job.runtime.backend in {"apptainer", "slurm"}:
        checks.append(
            Check(
                "Apptainer",
                shutil.which("apptainer") is not None,
                shutil.which("apptainer") or "not found",
            )
        )
        assert job.runtime.image is not None
        checks.append(_path_check("Container image", job.runtime.image, "Apptainer image"))
    if job.runtime.backend == "slurm":
        checks.append(
            Check("Slurm", shutil.which("sbatch") is not None, shutil.which("sbatch") or "not found")
        )
    return checks


def checkpoint_has_nan(path: Path) -> bool | None:
    """Return whether a NumPy checkpoint contains NaNs; unknown formats return None."""
    if path.suffix not in {".npy", ".npz"}:
        return None
    try:
        import numpy as np

        value = np.load(path)
        arrays = value.values() if hasattr(value, "values") else [value]
        return any(bool(np.isnan(array).any()) for array in arrays)
    except (ImportError, OSError, TypeError, ValueError):
        return None


def all_ok(checks: list[Check]) -> bool:
    return all(check.ok for check in checks)

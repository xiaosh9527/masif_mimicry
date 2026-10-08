from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class CheckpointReport:
    prefix: Path
    tensors: int
    floating_values: int
    nonfinite_values: int

    @property
    def ok(self) -> bool:
        return self.tensors > 0 and self.nonfinite_values == 0


def normalize_checkpoint_prefix(path: Path) -> Path:
    value = str(path.expanduser().resolve())
    for suffix in (".index", ".meta"):
        value = value.removesuffix(suffix)
    marker = ".data-"
    if marker in value:
        value = value.split(marker, 1)[0]
    return Path(value)


def inspect_checkpoint(path: Path) -> CheckpointReport:
    prefix = normalize_checkpoint_prefix(path)
    if not prefix.with_suffix(".index").exists():
        raise FileNotFoundError(f"TensorFlow checkpoint index not found: {prefix}.index")

    # TensorFlow 1.12 only exists in the native Python 3.6 environment. Keep
    # it out of the Python 3.11 CLI process and ask the legacy interpreter for
    # a small JSON report instead.
    masif_python = os.getenv("MASIF_PYTHON", shutil.which("python") or sys.executable)
    probe = """
import json
import sys
import numpy as np
import tensorflow as tf
reader = tf.train.load_checkpoint(sys.argv[1])
names = reader.get_variable_to_shape_map()
floating = 0
nonfinite = 0
for name in names:
    value = reader.get_tensor(name)
    if np.issubdtype(value.dtype, np.floating):
        floating += int(value.size)
        nonfinite += int((~np.isfinite(value)).sum())
print(json.dumps({'tensors': len(names), 'floating': floating, 'nonfinite': nonfinite}))
"""
    env = dict(os.environ)
    env["TF_CPP_MIN_LOG_LEVEL"] = "3"
    try:
        result = subprocess.run(
            [masif_python, "-c", probe, str(prefix)],
            check=True,
            capture_output=True,
            text=True,
            env=env,
        )
        values = json.loads(result.stdout.strip().splitlines()[-1])
    except (OSError, subprocess.CalledProcessError, json.JSONDecodeError, IndexError) as exc:
        detail = exc.stderr.strip() if isinstance(exc, subprocess.CalledProcessError) else str(exc)
        raise RuntimeError(
            f"Checkpoint inspection failed with MaSIF Python ({masif_python}): {detail}"
        ) from exc
    return CheckpointReport(
        prefix=prefix,
        tensors=int(values["tensors"]),
        floating_values=int(values["floating"]),
        nonfinite_values=int(values["nonfinite"]),
    )

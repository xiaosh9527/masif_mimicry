from pathlib import Path

import pytest

from masif_mimicry.checkpoint import normalize_checkpoint_prefix


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("model", "model"),
        ("model.index", "model"),
        ("model.meta", "model"),
        ("model.data-00000-of-00001", "model"),
    ],
)
def test_normalize_checkpoint_prefix(tmp_path: Path, value: str, expected: str):
    assert normalize_checkpoint_prefix(tmp_path / value) == tmp_path / expected

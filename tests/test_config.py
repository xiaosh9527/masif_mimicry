from pathlib import Path

import pytest
from pydantic import ValidationError

from masif_mimicry.config import JobConfig, Search, SearchPreset


@pytest.mark.parametrize(
    ("preset", "downsample"),
    [(SearchPreset.fast, 5), (SearchPreset.balanced, 3), (SearchPreset.sensitive, 1)],
)
def test_search_presets(preset: SearchPreset, downsample: int):
    assert Search(preset=preset).resolved()["downsample"] == downsample


def test_search_override_wins():
    assert Search(preset="fast", downsample=2).resolved()["downsample"] == 2


def test_target_processed_id():
    job = JobConfig.model_validate(
        {
            "name": "test",
            "target": {
                "structure": "/workspace/target.pdb",
                "pdb_id": "1abc",
                "query_chain": "A",
                "receptor_chain": "B",
            },
            "database": {"proteins": ["seed_A"]},
            "runtime": {"backend": "local", "suite_root": "/workspace"},
        }
    )
    assert job.target.processed_id == "1abc_A_B"


def test_container_paths_must_be_below_suite_root(tmp_path: Path):
    with pytest.raises(ValidationError, match="must be below"):
        JobConfig.model_validate(
            {
                "name": "test",
                "target": {"structure": "/outside/target.pdb", "query_chain": "A"},
                "database": {"proteins": ["seed_A"]},
                "runtime": {
                    "backend": "apptainer",
                    "suite_root": str(tmp_path),
                    "image": str(tmp_path / "image.sif"),
                },
            }
        )


def test_symlinked_manifest_is_allowed_when_real_root_is_bound(tmp_path: Path):
    suite = tmp_path / "scratch" / "masif_seed"
    database = tmp_path / "work" / "human_proteome"
    database.mkdir(parents=True)
    manifest = database / "hp_list.txt"
    manifest.write_text("seed_A\n", encoding="utf-8")
    link = suite / "masif" / "data" / "masif_human_proteome"
    link.parent.mkdir(parents=True)
    link.symlink_to(database, target_is_directory=True)

    job = JobConfig.model_validate(
        {
            "name": "array",
            "target": {
                "structure": str(suite / "masif_mimicry" / "target.pdb"),
                "query_chain": "A",
            },
            "database": {"protein_list": str(link / "hp_list.txt")},
            "postprocessing": {"enabled": False},
            "runtime": {
                "backend": "slurm",
                "suite_root": str(suite),
                "image": str(suite / "image.sif"),
                "binds": [str(database)],
                "slurm": {"array_tasks": 1},
            },
        }
    )

    assert job.database.protein_list == link / "hp_list.txt"


def test_database_requires_exactly_one_protein_source():
    base = {
        "name": "test",
        "target": {"structure": "/workspace/target.pdb", "query_chain": "A"},
        "runtime": {"backend": "local", "suite_root": "/workspace"},
    }
    with pytest.raises(ValidationError, match="exactly one"):
        JobConfig.model_validate({**base, "database": {}})
    with pytest.raises(ValidationError, match="exactly one"):
        JobConfig.model_validate(
            {
                **base,
                "database": {
                    "proteins": ["seed_A"],
                    "protein_list": "/workspace/seeds.txt",
                },
            }
        )


def test_slurm_arrays_require_search_only_mode(tmp_path: Path):
    suite = tmp_path / "suite"
    with pytest.raises(ValidationError, match="postprocessing must be disabled"):
        JobConfig.model_validate(
            {
                "name": "array",
                "target": {
                    "structure": str(suite / "masif_mimicry" / "target.pdb"),
                    "query_chain": "A",
                },
                "database": {"proteins": ["seed_A"]},
                "runtime": {
                    "backend": "slurm",
                    "suite_root": str(suite),
                    "image": str(suite / "image.sif"),
                    "slurm": {"array_tasks": 1},
                },
            }
        )

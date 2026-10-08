from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, Field, model_validator


class SearchPreset(str, Enum):
    fast = "fast"
    balanced = "balanced"
    sensitive = "sensitive"


PRESETS: dict[SearchPreset, dict[str, float | int]] = {
    SearchPreset.fast: {
        "downsample": 5,
        "descriptor_distance_cutoff": 1.5,
        "mimicry_score_cutoff": 0.45,
    },
    SearchPreset.balanced: {
        "downsample": 3,
        "descriptor_distance_cutoff": 2.0,
        "mimicry_score_cutoff": 0.4,
    },
    SearchPreset.sensitive: {
        "downsample": 1,
        "descriptor_distance_cutoff": 2.0,
        "mimicry_score_cutoff": 0.35,
    },
}


class Focus(BaseModel):
    residue: int = Field(gt=0)
    atom: str = Field(default="CA", min_length=1)


class Target(BaseModel):
    structure: Path
    pdb_id: str | None = None
    query_chain: str = Field(min_length=1)
    receptor_chain: str | None = Field(default=None, min_length=1)
    focus: Focus | None = None
    ligand: str | None = Field(default=None, pattern=r"^[^:]+:[^:]+$")

    @model_validator(mode="after")
    def infer_pdb_id(self) -> Target:
        if self.pdb_id is None:
            self.pdb_id = self.structure.stem
        return self

    @property
    def processed_id(self) -> str:
        chains = self.query_chain
        if self.receptor_chain:
            chains += f"_{self.receptor_chain}"
        return f"{self.pdb_id}_{chains}"


class Database(BaseModel):
    name: str = "masif_human_proteome"
    proteins: list[str] = Field(min_length=1)
    root: Path | None = None


class Search(BaseModel):
    preset: SearchPreset = SearchPreset.balanced
    downsample: int | None = Field(default=None, gt=0)
    descriptor_distance_cutoff: float | None = Field(default=None, gt=0)
    mimicry_score_cutoff: float | None = Field(default=None, ge=0, le=1)
    interface_cutoff: float = Field(default=0.0, ge=0, le=1)
    top_interface_percent: float = Field(default=0.0, ge=0, le=100)
    number_of_points: int = Field(default=5, gt=0)
    interface_only: bool = False
    count_clashes: bool = True
    ca_clash_threshold: float = Field(default=1.0, ge=0)
    heavy_atom_clash_threshold: float = Field(default=5.0, ge=0)

    def resolved(self) -> dict[str, Any]:
        values = dict(PRESETS[self.preset])
        if self.downsample is not None:
            values["downsample"] = self.downsample
        if self.descriptor_distance_cutoff is not None:
            values["descriptor_distance_cutoff"] = self.descriptor_distance_cutoff
        if self.mimicry_score_cutoff is not None:
            values["mimicry_score_cutoff"] = self.mimicry_score_cutoff
        values.update(
            interface_cutoff=self.interface_cutoff,
            top_interface_percent=self.top_interface_percent,
            number_of_points=self.number_of_points,
            interface_only=self.interface_only,
            count_clashes=self.count_clashes,
            ca_clash_threshold=self.ca_clash_threshold,
            heavy_atom_clash_threshold=self.heavy_atom_clash_threshold,
        )
        return values


class Postprocessing(BaseModel):
    enabled: bool = True
    truncation_length: int = Field(default=86, gt=0)


class SlurmResources(BaseModel):
    partition: str | None = None
    cpus: int = Field(default=4, gt=0)
    memory: str = "24G"
    time: str = Field(default="24:00:00", pattern=r"^\d{1,3}:\d{2}:\d{2}$")


class Runtime(BaseModel):
    backend: Literal["local", "apptainer", "slurm"] = "slurm"
    suite_root: Path
    image: Path | None = None
    submit: bool = False
    slurm: SlurmResources = Field(default_factory=SlurmResources)

    @model_validator(mode="after")
    def require_image_for_container_backends(self) -> Runtime:
        if self.backend in {"apptainer", "slurm"} and self.image is None:
            raise ValueError(f"runtime.image is required for the {self.backend} backend")
        return self


class JobConfig(BaseModel):
    name: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
    target: Target
    database: Database
    search: Search = Field(default_factory=Search)
    postprocessing: Postprocessing = Field(default_factory=Postprocessing)
    runtime: Runtime
    output_root: Path = Path("runs")

    @model_validator(mode="after")
    def validate_container_paths(self) -> JobConfig:
        if self.runtime.backend == "local":
            return self
        suite = self.runtime.suite_root.resolve()
        paths = {
            "target.structure": self.target.structure,
            "output_root": self.output_root,
        }
        if self.database.root is not None:
            paths["database.root"] = self.database.root
        for label, path in paths.items():
            if path.is_absolute() and not path.resolve().is_relative_to(suite):
                raise ValueError(
                    f"{label} ({path}) must be below runtime.suite_root ({suite}) "
                    "for container backends"
                )
        return self

    def resolve_paths(self, base_dir: Path) -> JobConfig:
        data = self.model_dump()
        for section, key in (
            ("target", "structure"),
            ("database", "root"),
            ("runtime", "suite_root"),
            ("runtime", "image"),
        ):
            value = data[section].get(key)
            if value is not None:
                path = Path(value).expanduser()
                data[section][key] = path if path.is_absolute() else (base_dir / path).resolve()
        output = Path(data["output_root"]).expanduser()
        data["output_root"] = output if output.is_absolute() else (base_dir / output).resolve()
        return JobConfig.model_validate(data)

    @property
    def run_dir(self) -> Path:
        return self.output_root / self.name


def load_job(path: Path) -> JobConfig:
    path = path.expanduser().resolve()
    with path.open(encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    if not isinstance(raw, dict):
        raise TypeError("Job configuration must be a YAML mapping")
    return JobConfig.model_validate(raw).resolve_paths(path.parent)


def dump_job(job: JobConfig, path: Path) -> None:
    data = job.model_dump(mode="json", exclude_none=True)
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")

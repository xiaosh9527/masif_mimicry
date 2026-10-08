# MaSIF-mimicry

MaSIF-mimicry searches protein surfaces for patches that resemble a target interface. This repository contains the workflow used for the experiments in [*Proteome-wide identification of the druggable CRBN interactome*](https://www.nature.com/articles/s41587-026-03237-7).

![MaSIF-mimicry overview and pipeline](img/masif_mimicry.jpg)

> [!IMPORTANT]
> MaSIF-mimicry currently depends on the MaSIF and MaSIF-seed research environments and has been tested on x86 Linux. Run `masif-mimicry doctor job.yaml` before starting a long job. Also verify that your downloaded neural-network checkpoint is finite; see the open [checkpoint reproducibility report](https://github.com/xiaosh9527/masif_mimicry/issues/3).

## Quick start

Clone this repository beside the two upstream projects:

```text
masif_seed/
├── masif/
├── masif_seed_search/
└── masif_mimicry/
```

Build the single runtime image:

```bash
cd masif_mimicry
docker build --platform linux/amd64 -t masif-mimicry .
```

The image contains two deliberately isolated interpreters—native Python 3.6
for MaSIF and Conda Python 3.11 for bonsai and the CLI. No host Python
environment or post-build Conda installation is required.

For Apptainer, build from the Docker image or publish the image and pull it:

```bash
apptainer build masif_mimicry.sif docker-daemon://masif-mimicry:latest
```

Create an experiment:

```bash
apptainer exec masif_mimicry.sif masif-mimicry init my-screen
cd my-screen
```

Put the target PDB in `inputs/`, edit the generated `job.yaml`, and validate the setup:

```bash
apptainer exec --bind /path/to/masif_seed:/workspace masif_mimicry.sif \
  masif-mimicry validate job.yaml
apptainer exec --bind /path/to/masif_seed:/workspace masif_mimicry.sif \
  masif-mimicry doctor job.yaml
apptainer exec --bind /path/to/masif_seed:/workspace masif_mimicry.sif \
  masif-mimicry plan job.yaml
```

Validate downloaded TensorFlow weights from inside the MaSIF environment before a large run:

```bash
apptainer exec --bind /path/to/masif_seed:/workspace masif_mimicry.sif \
  masif-mimicry checkpoint /workspace/masif/data/.../model
```

Generate and run the workflow:

```bash
apptainer exec --bind /path/to/masif_seed:/workspace masif_mimicry.sif \
  masif-mimicry run job.yaml
```

For Slurm, the default is intentionally safe: it creates and validates `submit.slurm` without submitting it. Review the file, then either run `sbatch` yourself or set `runtime.submit: true`.

Inspect progress and create a browser-readable report:

```bash
masif-mimicry status runs/znf692-vs-mtor
masif-mimicry results runs/znf692-vs-mtor --open
```

The form-based browser interface is also installed in the image:

```bash
apptainer exec --bind /path/to/masif_seed:/workspace masif_mimicry.sif \
  masif-mimicry ui
```

For development outside the image, the lightweight CLI can still be installed
on the host. This is optional and is not part of the production runtime:

```bash
python -m pip install -e '.[dev]'
```

## One configuration file

The interface replaces editable shell scripts with a validated YAML job:

```yaml
name: znf692-vs-mtor

target:
  structure: inputs/6h0g.pdb
  pdb_id: 6h0g
  query_chain: C
  receptor_chain: B
  focus: {residue: 423, atom: CA}
  ligand: B:Y70

database:
  name: masif_human_proteome
  proteins: [P42345-F1-dom-01_A]

search:
  preset: balanced
  count_clashes: true

postprocessing:
  enabled: true
  truncation_length: 86

runtime:
  backend: local
  suite_root: /workspace
  submit: false
  slurm:
    cpus: 4
    memory: 24G
    time: "24:00:00"

output_root: runs
```

Paths are resolved relative to the YAML file. When the CLI runs inside the
container, use `backend: local` and `/workspace` as `suite_root`. The local
backend means “already inside the image”; it does not combine the two Python
environments.

## Python compatibility

The image has exactly two Python runtimes:

| Path | Version | Used for |
| --- | --- | --- |
| `/usr/local/bin/masif-python` | 3.6 | Surface preparation, TensorFlow 1.12, MaSIF search |
| `/opt/conda/envs/bonsai/bin/python` | 3.11 | CLI, RDKit, postprocessing, EvoEF2 orchestration |

Generated workflows call these absolute paths. They never run `conda activate`,
and Miniconda is deliberately not added to the global `PATH`. Environment
variables `MASIF_PYTHON` and `BONSAI_PYTHON` provide equivalent explicit paths
for the legacy scripts.

## Search presets

| Preset | Downsample | Descriptor cutoff | Mimicry-score cutoff | Intended use |
| --- | ---: | ---: | ---: | --- |
| `fast` | 5 | 1.5 | 0.45 | Quick exploratory screen |
| `balanced` | 3 | 2.0 | 0.40 | Recommended starting point |
| `sensitive` | 1 | 2.0 | 0.35 | More candidates and more computation |

Individual values can override a preset:

```yaml
search:
  preset: balanced
  downsample: 2
  descriptor_distance_cutoff: 1.8
  mimicry_score_cutoff: 0.42
```

## Runtime backends

- `local` runs directly in the current environment and is the normal choice when the CLI is already running inside the image.
- `apptainer` runs the workflow in the configured image and mounts the suite at `/workspace`.
- `slurm` generates a reproducible Slurm submission script which runs the Apptainer backend.

Large databases are expected to be preprocessed. Each value in `database.proteins` is a processed MaSIF identifier such as `P42345-F1-dom-01_A`.

Target preprocessing is cached globally under
`masif_mimicry/data/target/targets/<processed-id>`. A later run using the same
processed target ID reuses its surface features and descriptors. Each run only
contains a symlink at `data/target`, keeping older status and inspection tools
compatible without copying the target data. To deliberately rebuild a target,
delete its directory from the shared cache before starting the run.

## Run directory

Every experiment has a predictable structure:

```text
runs/<name>/
├── run.yaml             # fully resolved configuration
├── provenance.json      # interface, host, time, and backend
├── run.sh|submit.slurm  # generated execution script
├── inputs/
├── data/target -> shared target cache
├── results/
├── postprocess/
├── truncate-86/
├── logs/
└── report.html
```

The generated `run.yaml`, execution script, and provenance file should be archived with published results.

## Existing installations

The original scripts remain available for compatibility:

- `preprocess_pdb.sh` preprocesses an individual PDB.
- `source/masif_mimicry_search.py` runs the core search.
- `scripts/bonsai_scripts/` contains postprocessing and truncation algorithms.
- `data/template/run.slurm` is the original worked example.

New experiments should use the `masif-mimicry` command so parameters, paths, and outputs are recorded consistently.

## Development

```bash
python -m pip install -e '.[dev]'
pytest
ruff check src tests
```

The interface tests do not run the multi-gigabyte scientific pipeline. They validate configuration, parameter presets, path translation, generated execution scripts, and run metadata.

## Hardware and dependencies

Surface preprocessing invokes external tools including MSMS, APBS, PDB2PQR,
Reduce, and USalign. Search uses the image's native Python 3.6 environment.
Bonsai uses the embedded Conda Python 3.11 environment together with EvoEF2,
STRIDE, and the bundled DeepTMHMM predictions.

APBS surface-charge calculation is the main preprocessing bottleneck. A distributed CPU cluster is recommended for proteome-scale preprocessing; a single prepared protein can be used for the included example.

## License and citation

MaSIF-mimicry is released under the Apache License 2.0. If you use it, cite [*Proteome-wide identification of the druggable CRBN interactome*](https://www.nature.com/articles/s41587-026-03237-7).

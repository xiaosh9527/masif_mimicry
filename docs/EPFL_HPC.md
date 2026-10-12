# EPFL/LPDI HPC notes

> This document is for EPFL/LPDI users. It contains site-specific filesystem
> paths and shell settings that should not be copied to other clusters.

## Storage layout

The working checkout is:

```text
/scratch/shxiao/masif_seed
```

The human-proteome path inside that checkout is a symbolic link whose real
location is:

```text
/work/lpdi/users/shxiao/masif_seed/masif/data/masif_human_proteome_domains_merged
```

The real location must be listed under `runtime.binds`; the launcher then
mounts it automatically in every array task.

## One-time shell setup

```bash
export MASIF_SEED_ROOT="/scratch/shxiao/masif_seed"
export MIMICRY_IMAGE="$MASIF_SEED_ROOT/masif_mimicry/masif_mimicry.sif"
export PATH="$MASIF_SEED_ROOT/masif_mimicry/bin:$PATH"
```

These lines may be added to `~/.bashrc` for future sessions.

## EPFL array configuration

Save the following as
`$MASIF_SEED_ROOT/masif_mimicry/job.epfl.yaml` in the cluster checkout:

```yaml
name: whole-human-proteome

target:
  structure: /scratch/shxiao/masif_seed/masif_mimicry/data/template/input/6h0g.pdb
  pdb_id: 6h0g
  query_chain: C
  receptor_chain: B
  focus: {residue: 423, atom: CA}
  ligand: B:Y70

database:
  name: masif_human_proteome
  protein_list: /scratch/shxiao/masif_seed/masif/data/masif_human_proteome/lists/hp_list.txt

search:
  preset: balanced
  count_clashes: true

postprocessing:
  enabled: false

runtime:
  backend: slurm
  suite_root: /scratch/shxiao/masif_seed
  image: /scratch/shxiao/masif_seed/masif_mimicry/masif_mimicry.sif
  binds:
    - /work/lpdi/users/shxiao/masif_seed/masif/data/masif_human_proteome_domains_merged
  submit: true
  slurm:
    nodes: 1
    qos: serial
    ntasks_per_node: 1
    job_name: benchmark
    array_tasks: 5
    cpus: 2
    memory: 10g
    time: "24:00:00"
    output: logs/%x-%A_%a.out
    # array_max_parallel: 50

output_root: /scratch/shxiao/masif_seed/masif_mimicry/runs
```

Use five tasks for the initial cluster test. Increase `array_tasks` only after
the test succeeds and according to the current Slurm policy.

## Submit, monitor, and resume

No manual Apptainer or `sbatch` command is needed:

```bash
masif-mimicry run \
  "$MASIF_SEED_ROOT/masif_mimicry/job.epfl.yaml"
```

Check progress:

```bash
masif-mimicry status \
  "$MASIF_SEED_ROOT/masif_mimicry/runs/whole-human-proteome"
```

Repeat the same `run` command to resume. Completed shards are skipped.

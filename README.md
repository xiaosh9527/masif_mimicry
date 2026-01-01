## _MaSIF-mimicry_ 

This repository contains code to perform the experiments in _Mapping the latent CRBN-molecular glue degrader interactome_

### Pipeline overview
![MaSIF-mimicry overview and pipeline](img/masif_mimicry.jpg)

## Table of Contents: 
- [System requirements](#system-requirements)
  - [Hardware requirements](#hardware-requirements)
  - [Software requirements](#software-requirements)
- [Installation](#installation)
- [Step-by-step example with Docker](#step-by-step-example-with-docker)
- [Configuring parameters](#configuring-parameters)
- [License](#license)
- [Reference](#reference)

## System requirements

### Hardware requirements

MaSIF-seed has been tested on Linux, and it is recommended to run on an x86-based Linux Docker container. 

Currently, MaSIF takes a few seconds to preprocess every protein. We find the main bottleneck is the APBS computation of surface charges, which can likely be optimized. Nevertheless, we recommend a distributed CPU cluster for preprocessing large protein datasets.

### Software requirements

MaSIF relies on external software/libraries to handle protein databank files and surface files, 
to compute chemical/geometric features and coordinates, and to perform neural network calculations. 
These dependencies have been installed in an Apptainer image. 
The pipeline also requires additional packages for performing the postprocessing _bonsai_ workflow. 
An environment.yml file is provided to create a conda environment that is independent from the MaSIF env.
Below, we give an example of how to install these separately in an Apptainer image. 

## Installation

MaSIF-mimicry has been tested on Linux. To run the mimicry pipeline, first clone the official MaSIF-mimicry repository and then clone this repository inside it. 

MaSIF is written in Python and does not require compilation. Since MaSIF relies on a few external programs (MSMS, APBS) and libraries (PyMesh, TensorFlow, SciPy, Open3D), we strongly recommend that you use the Dockerfile and Docker container. Setting up the environment should take only a few minutes.
```
git clone https://github.com/LPDI-EPFL/masif_seed.git
cd masif_seed
git clone https://github.com/xiaosh9527/masif_mimicry.git
cd masif_mimicry
```

We use Apptainer/Singularity on an HPC system. Apptainer image is created as follows:

```
apptainer pull masif_mimicry.sif docker://shxiao/masif_mimicry:latest
```

## Step-by-step example
We provide an example of searching for surface patches in P42345 (mTOR) that mimic the 6h0g_C (ZNF692) degron interface.
We first divide proteins in DB (i.e., P42345) into structured domains, and use MaSIF to process these domains as well as the target protein, chain C of 6h0g ZNF692. Explanation for each step is detailed in the Slurm file. 
```
sbatch preprocess_pdb_apptainer.slurm
```

Once both the DB and target protein are processed, run the script to search for a surface patch in mTOR that mimics the ZNF692 degron interface. 

```
cd ./data/template
sbatch run.slurm
```

Before proceeding to the postprocessing step (_bonsai_), we create a conda env for the required dependencies (EvoEF2, Stride, DeepTMHMM) in the Apptainer image. 
This is to ensure all dependencies are correctly installed and configured as a separate environment from the MaSIF environment.
```
cd ../../scripts/bonsai_scripts/
bash ./slurm/conda.sh
```

Ideally, one should collect all hits from the DB and create a directory that contains PDBs of the best poses. Here, we use the raw pose output from the P42345 (mTOR) vs 6h0g_C (ZNF692) search. 
Finally, we run the postprocessing script to truncate the hits and compute metrics. 

```
bash postprocessing_truncation_workflow_apptainer.sh
```

By default, the output metrics are stored in the `Truncate_86/proc_truncated_scores.csv` file inside the defined output directory, and all the truncated PDBs are stored in the `Truncate_86/truncate` directory.

## Configuring parameters

Some parameters may improve your search:

The main criterion for speed is the 'descriptor distance cutoff'. This cutoff determines which fingerprints are further considered and which are completely discarded. The lower the value, the faster the search (and the higher the number of false negatives):
```
# Score cutoffs -- these are empirical values; if they are too loose, then you get a lot of results.
# Descriptor distance cutoff for the patch. All scores below this value are accepted for further processing.
desc_dist_cutoff = 1.5 # Recommended values: [1.5-2.0] (lower is stricter)
```

Another significant cutoff value is the interface cutoff. All patches are scored due to their interface propensity. One can assume that during the search, we want to find peptides with a high interface cutoff. You can lower this value to increase the number of candidates:

```
# Interface cutoff value: all patches with a score below this value are discarded. For the degron interface, as MaSIF has not been trained with ternary complexes, we recommend an iface cutoff at 0.0
iface_cutoff = 0.0 # Recommended values: [0.0-0.95] range (higher is stricter)
```

Finally, for the normalized desc_dist_score (i.e., MaSIF mimicry score), we recommend a cut-off at 0.4:

```
# MaSIF mimicry score cutoff - Discard anything below this score
desc_dist_score_cutoff = 0.4 # Recommended values: [0.35-0.5] (higher is stricter)
```

The number of sites to target in the protein (generally 1 should be fine):

```
# Number of sites to target in the protein
num_sites = 10
```

Due to the large size of the proteome dataset, it is possible to perform searches with downsampling and select every k patches:

```
# Downsampling every 3 patches
downsample = 3
```

## License

MaSIF-mimicry is released under an [Apache v2.0 license](LICENSE).

## Reference
If you use this code, please use the bibtex entry in [citation.bib](citation.bib)

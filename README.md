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

MaSIF-seed has been tested on Linux, and it is recommended to run on an x86-based linux Docker container. 

Currently, MaSIF takes a few seconds to preprocess every protein. We find the main bottleneck to be the APBS computation for surface charges, which can likely be optimized. Nevertheless, we recommend a distributed CPU cluster to preprocess the data for large datasets of proteins.

### Software requirements

MaSIF relies on external software/libraries to handle protein databank files and surface files, 
to compute chemical/geometric features and coordinates, and to perform neural network calculations. 
The following is the list of required libraries and programs, as well as the version on which it was tested (in parentheses).
* [Python](https://www.python.org/) (3.6)
* [reduce](http://kinemage.biochem.duke.edu/software/reduce.php) (3.23). To add protons to proteins. 
* [MSMS](http://mgltools.scripps.edu/packages/MSMS/) (2.6.1). To compute the surface of proteins. 
* [BioPython](https://github.com/biopython/biopython) (1.66). To parse PDB files. 
* [PyMesh](https://github.com/PyMesh/PyMesh) (0.1.14). To handle ply surface files, attributes, and to regularize meshes.
* PDB2PQR (2.1.1), multivalue, and [APBS](http://www.poissonboltzmann.org/) (1.5). These programs are necessary to compute electrostatics charges.
* [Open3D](https://github.com/IntelVCL/Open3D) (0.5.0.0). Mainly used for RANSAC alignment.
* [Tensorflow](https://www.tensorflow.org/) (1.9). Use to model, train, and evaluate the actual neural networks. Models were trained and evaluated on a NVIDIA Tesla K40 GPU.
* [Pymol](https://pymol.org/2/) (2.5.0). This optional program allows one to visualize surface files.
* [python-igraph](https://igraph.org/python/) (0.9.6). Used in pae_to_domain package to split AlphaFold (AF) models into individual domains. Only needed if you want to process AF models.
* [pae_to_domain](https://github.com/tristanic/pae_to_domains). Code used to split AF models into individual domains based on Predicted Aligned Error (PAE) values. Only needed if you want to process AF models.

## Installation

MaSIF-mimicry has been tested on Linux. To run the mimicry pipeline, first clone the official MaSIF-mimicry repository and then clone this repository inside it. 

MaSIF is written in Python and does not require compilation. Since MaSIF relies on a few external programs (MSMS, APBS) and libraries (PyMesh, Tensorflow, Scipy, Open3D), we strongly recommend you use the Dockerfile and Docker container. Setting up the environment should take a few minutes only.
```bash
git clone https://github.com/LPDI-EPFL/masif_seed.git
cd masif_seed
git clone https://github.com/xiaosh9527/masif_mimicry.git
cd masif_mimicry
```

We use Apptainer/Singularity on an HPC system. 
Apptainer image is created as follows:

```
apptainer pull masif_mimicry.sif docker://shxiao/masif_mimicry:latest
```

## Step-by-step example
We provide an example of searching for surface patches in P42345 (mTOR) that mimic the 6h0g_C (ZNF692) degron interface.
Explaination for each step is detailed in the slurm file. 
```
sbatch preprocess_pdb_apptainer.slurm
```

Once both the DB and target protein are processed, run the script to search for a surface patch in mTOR that mimics the ZNF692 degron interface. 

```
cd ./data/template
sbatch run.slurm
```

Before proceed to the postprocessing step (_bonsai_), we create a conda env for the required dependencies (EvoEF2, Stride, DeepTMHMM) in the Apptainer image. 
This is to ensure all dependencies are properly installed and configured as a separate environment from the MaSIF environment.
```bash
cd ../../scripts/bonsai_scripts/
bash ./slurm/conda.sh
```

Ideally one should collect all hits from the DB and create a directory that contains PDBs of the best poses. Here we use the raw output of poses from P42345 (mTOR) vs 6h0g_C (ZNF692) search. 
Finally, we run the postprocessing script to truncate the hits and compute metrics. 

```bash
bash postprocessing_truncation_workflow_apptainer.sh
```

By default, the output metrics are stored in `Truncate_86/proc_truncated_scores.csv` file inside of the defined output directory, and all the truncated PDBs are stored in the `Truncate_86/truncate` directory.

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



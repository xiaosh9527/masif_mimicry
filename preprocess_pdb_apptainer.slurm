#!/bin/bash -l
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 24G
#SBATCH --time 00:30:00
#SBATCH --output logs/masif_%J.out

IMAGE=./masif_mimicry.sif
MASIF_DIR=../

##########################################################
# Preprocess human proteome structures into domains
##########################################################

UNIPROT_ID="P42345"  # Example: mTOR
MASIF_DB_NAME="masif_human_proteome"  # Directory containing MaSIF human proteome data

# To reproduce the experiments in the paper, the entire datasets for the human proteome consume several terabytes. We will test MaSIF-mimicry using one example between mTOR and Ikaros.
# Here, we provide an example to process mTOR structures into individual structural domains by using pae_to_domain to split the AlphaFold models based on a PAE cutoff of 15 Å.

srun apptainer exec --bind $MASIF_DIR:/workspace $IMAGE /bin/bash -c \
       "cd /workspace/masif/data/ &&
       mkdir -p $MASIF_DB_NAME &&
       cd $MASIF_DB_NAME &&
       if [ ! -d "pae_to_domains" ] ; then
            git clone https://github.com/tristanic/pae_to_domains.git
       fi &&
       sh ../../../masif_mimicry/scripts/process_af_model.sh $UNIPROT_ID"

# This will generate truncated structural domains from the original full-length AF model. The same process can be applied to rest of the proteome proteins to generate a MaSIF human proteome database. 
# Note: For the article we used AF database v4, which can be downladed in this [link](https://ftp.ebi.ac.uk/pub/databases/alphafold/v4/UP000005640_9606_HUMAN_v4.tar) (4.8 GB).

# Warning: The current script [`domain_split.py`](scripts/domain_split.py), called by [`process_af_model.sh`](scripts/process_af_model.sh), uses the version v6 because v4 is no longer supported for download using the Uniprot ID.
# If you'd like to preprocess more than one protein, please modify the script with array jobs or loops.

##############################################
# Preprocess AF domain DB
##############################################
# Now that we have the structural domains, we can preprocess them using the MaSIF pipeline.

# Below, we use one of the mTOR domains as an example. To process a fragment we will use the [`preprocess_pdb.sh`](preprocess_pdb.sh) script, which runs the following tasks:
# 1. Preprocess the PDB file (add hydrogens, compute surface with MSMS, compute geometric and chemical features...)
# 2. Precompute the MaSIF-site features
# 3. Precompute the MaSIF-search features
# 4. Predict interface sites with MaSIF-site
# 5. Compute MaSIF descriptors
# 6. Copy files required for mimicry search to output directory

srun apptainer exec --bind $MASIF_DIR:/workspace $IMAGE /bin/bash -c \
       "cd /workspace/masif/data/masif_human_proteome &&
       bash /workspace/masif_mimicry/preprocess_pdb.sh input/fragments/$UNIPROT_ID-F1-dom-01.pdb $UNIPROT_ID-F1-dom-01_A -o /workspace/masif/data/masif_human_proteome"

# Note: In cases where multiple domains are generated, the same procedure must be applied to all.

################################################
# Preprocess target PDB
################################################

# Once the site predictions and descriptors on `P42345-F1-dom-01_A` have been computed, we can focus on the target. 
# We use the provided example target structure `6h0g.pdb` (ZNF692 in complex with CRBN and ligand B_Y70).

PDB_ID="input/6h0g.pdb"  # Example: ZNF692
PDB_CHAIN="6h0g_C_B"  # Chain C (ZNF692) and chain B (CRBN)
TARGET_OUTPUT_DIR="/workspace/masif_mimicry/data/template/"  # Output directory for target PDB

srun apptainer exec --bind $MASIF_DIR:/workspace $IMAGE /bin/bash -c \
        "cd $TARGET_OUTPUT_DIR &&
        bash /workspace/masif_mimicry/preprocess_pdb.sh $PDB_ID $PDB_CHAIN -o $TARGET_OUTPUT_DIR"
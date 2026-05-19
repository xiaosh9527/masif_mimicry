#!/bin/bash -l
#SBATCH --job-name=preprocess
#SBATCH --partition=standard
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 7000
#SBATCH --time 00:30:00
#SBATCH --output logs/masif_%J.out

IMAGE=./masif_mimicry.sif
MASIF_DIR=../

################################################
# Preprocess target PDB
################################################

# Once the site predictions and descriptors on `P42345-F1-dom-01_A` have been computed, we can focus on the target. 
# We use the provided example target structure `6h0g.pdb` (ZNF692 in complex with CRBN and ligand B_Y70).

PDB_ID="input/021structure.pdb"  # Example: ZNF692
PDB_CHAIN="021structure_C_AB"  # p1: chain C; p2: chains A and B (partner structure)
TARGET_OUTPUT_DIR="/workspace/masif_mimicry/data/NUP98/"  # Output directory for target PDB
LIGAND="A_021"

srun apptainer exec --bind $MASIF_DIR:/workspace $IMAGE /bin/bash -c \
        "cd $TARGET_OUTPUT_DIR &&
        bash /workspace/masif_mimicry/preprocess_pdb.sh $PDB_ID $PDB_CHAIN -o $TARGET_OUTPUT_DIR -l $LIGAND -s input/$LIGAND.sdf"
#!/bin/bash -l
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 24G
#SBATCH --time 00:30:00
#SBATCH --output logs/masif_%J.out

IMAGE=./masif_mimicry.sif
MASIF_DIR=../

srun apptainer exec --bind $MASIF_DIR:/workspace $IMAGE /bin/bash -c \
       "cd /workspace/masif/data/ &&
       mkdir -p masif_human_proteome &&
       cd masif_human_proteome &&
       if [ ! -d "pae_to_domains" ] ; then
            git clone https://github.com/tristanic/pae_to_domains.git
       fi &&
       sh ../../../masif_mimicry/scripts/process_af_model.sh P42345"

srun apptainer exec --bind $MASIF_DIR:/workspace $IMAGE /bin/bash -c \
       "cd /workspace/masif/data/masif_human_proteome &&
       bash /workspace/masif_mimicry/preprocess_pdb.sh input/fragments/P42345-F1-dom-01.pdb P42345-F1-dom-01_A"

srun apptainer exec --bind $MASIF_DIR:/workspace $IMAGE /bin/bash -c \
        "cd /workspace/masif_mimicry/data/template/ &&
        bash /workspace/masif_mimicry/preprocess_pdb.sh input/6h0g.pdb 6h0g_C_B"
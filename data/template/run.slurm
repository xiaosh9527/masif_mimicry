#!/bin/bash -l
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 6G
#SBATCH --time 24:00:00
#SBATCH --output logs/run_%J.out

masif_root=/scratch/shxiao/masif_seed/
IMAGE=$masif_root/masif_mimicry/masif_mimicry.sif

apptainer exec \
  --bind "$masif_root:/workspace" \
  --pwd /workspace \
  "$IMAGE" /bin/bash -c '
    cd $masif_mimicry_root/data/template/
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8

    masif_db_root=/workspace/masif
    masif_seed_root=/workspace/masif_seed_search
    masif_mimicry_root=/workspace/masif_mimicry

    masif_db_source=$masif_db_root/source
    masif_seed_source=$masif_seed_root/source
    masif_mimicry_source=$masif_mimicry_root/source

    export PYTHONPATH="$PYTHONPATH:$masif_db_source:$masif_seed_source:$masif_mimicry_source:$PWD"

    python -u "$masif_mimicry_source/masif_mimicry_search.py" \
      --seed_db_root "$masif_db_root" \
      --seed_db masif_human_proteome \
      --target_db_root "$masif_mimicry_root" \
      --target_db template \
      --target_pdb 6h0g_C_B \
      --target_chain C \
      --target_residue 423 \
      --target_atom CA \
      --top_iface_percent 0.0 \
      --count_clashes \
      --downsample 5 \
      --desc_dist_cutoff 2.0 \
      --desc_dist_score_cutoff 0.4 \
      --output_dir "$PWD/results/" \
      --output_postfix test \
      --seed_pdb P42345-F1-dom-01_A
  '
#!/bin/bash

conda activate masif

masif_root=/work/lpdi/users/shxiao/masif_seed
masif_db_root=$masif_root/masif
masif_seed_root=$masif_root/masif_seed_search
masif_mimicry_root=$masif_root/masif_mimicry

masif_db_source=$masif_db_root/source
masif_seed_source=$masif_seed_root/source
masif_mimicry_source=$masif_mimicry_root/source

export masif_db_source
export masif_seed_source
export masif_mimicry_source
export PYTHONPATH=$PYTHONPATH:$masif_db_source:$masif_seed_source:$masif_mimicry_source:`pwd`

python -u $masif_mimicry_source/masif_mimicry_search.py \
--seed_db_root $masif_db_root \
--seed_db masif_human_proteome \
--target_db_root $masif_seed_root \
--target_db masif_mimicry_search \
--target_pdb 6h0g_C_B \
--target_chain C \
--target_residue 423 \
--target_atom CA \
--top_iface_percent 0.0 \
--count_clashes \
--downsample 5 \
--desc_dist_cutoff 2.0 \
--desc_dist_score_cutoff 0.4 \
--output_dir ./results/ \
--output_postfix test \
--seed_pdb P42345-F1-dom-01_A

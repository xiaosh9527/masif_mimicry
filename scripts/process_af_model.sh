#!/bin/bash

source /work/lpdi/users/shxiao/masif_seed/load_masif_environment_jed.sh

masif_root=/work/lpdi/users/shxiao/masif_seed
masif_mimicry_root=$masif_root/masif_mimicry
masif_mimicry_script=$masif_mimicry_root/scripts

export masif_mimicry_script
export PYTHONPATH=$PYTHONPATH:$masif_mimicry_script:`pwd`

python -u $masif_mimicry_script/domain_split.py \
--uniprot_id $1 \
--pae_cutoff 15.0 \
--input_dir ./input/
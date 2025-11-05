#!/bin/bash

# USER CHANGEABLE PATH
masif_root=/workspace # Assuming you're running the code in the docker container
# DO NOT CHANGE BELOW
masif_mimicry_root=$masif_root/masif_mimicry
masif_mimicry_script=$masif_mimicry_root/scripts

export masif_mimicry_script
export PYTHONPATH=$PYTHONPATH:$masif_mimicry_script:`pwd`

python -u $masif_mimicry_script/domain_split.py \
--uniprot_id $1 \
--pae_cutoff 15.0 \
--input_dir ./input/
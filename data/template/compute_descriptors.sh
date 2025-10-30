#!/bin/bash
# masif_root=$(git rev-parse --show-toplevel)
masif_root=/work/lpdi/users/shxiao/masif_seed
masif_source=$masif_root/masif/source/
masif_data=$masif_root/masif/data/
export PYTHONPATH=$PYTHONPATH:$masif_source:$masif_data/masif_ppi_search/
python $masif_source/masif_ppi_search/masif_ppi_search_comp_desc.py nn_models.sc05.all_feat.custom_params $1 $2

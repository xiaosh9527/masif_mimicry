#!/bin/bash -l
#SBATCH --job-name=mimicry_search
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 6G
#SBATCH --time 02:00:00
#SBATCH --array=0
#SBATCH --output logs/mimicry_%A/slurm_%A_%a.out
#SBATCH --error logs/mimicry_%A/slurm_%A_%a.out

set -euo pipefail

# User specific path:
masif_root=/scratch/ymeng/masif_seed

# Run specific paths:
DATABASE_DIR="/scratch/ymeng/TED_domainome/output"
DATABASE_SUBSET_DIR="/scratch/ymeng/TED_domainome/output/filtered_intracellular_domainome/subsets"

# DATABASE_DIR="/work/lpdi/users/shxiao/masif_seed/masif/data/masif_human_proteome_domains_merged"
# DATABASE_SUBSET_DIR="/work/lpdi/users/schneuing/tricomplex-design/search_lists/filtered20251016"
TARGET_PREPROCESS_DIR="data/NUP98"
TARGET_PDB="021structure_C_AB"
TARGET_CHAIN="C"
TARGET_RESIDUE=728
TARGET_ATOM="CB"

OUTPUT_DIR="data/NUP98/search_results"


# ------ Do not change below this line ------
masif_db_root="${masif_root}/masif"
masif_seed_root="${masif_root}/masif_seed_search"
masif_mimicry_root="${masif_root}/masif_mimicry"

IMAGE="${masif_mimicry_root}/masif_mimicry.sif"

masif_db_source="${masif_db_root}/source"
masif_seed_source="${masif_seed_root}/source"
masif_mimicry_source="${masif_mimicry_root}/source"

export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}${masif_db_source}:${masif_seed_source}:${masif_mimicry_source}"

apptainer exec \
  --bind /scratch:/scratch \
  --bind /work:/work \
  --env PYTHONPATH="${PYTHONPATH}" \
  "${IMAGE}" \
  python -u "${masif_mimicry_source}/masif_mimicry_search.py" \
    --database_dir "${DATABASE_DIR}" \
    --target_preprocess_dir "${TARGET_PREPROCESS_DIR}" \
    --target_pdb "${TARGET_PDB}" \
    --target_chain "${TARGET_CHAIN}" \
    --target_residue "${TARGET_RESIDUE}" \
    --target_atom "${TARGET_ATOM}" \
    --target_ppi_id "p1" \
    --top_iface_percent 0.0 \
    --count_clashes \
    --ca_clash_threshold 1.0 \
    --heavy_atom_clash_threshold 5.0 \
    --num_points 20 \
    --target_sampling_radius 5 \
    --desc_dist_cutoff 2.0 \
    --desc_dist_score_cutoff 0.45 \
    --output_dir "${OUTPUT_DIR}" \
    --split_seed_list ${DATABASE_SUBSET_DIR}/$SLURM_ARRAY_TASK_ID

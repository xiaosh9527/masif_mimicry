# Shared configuration for the MaSIF mimicry SLURM pipeline.
# Sourced by 1_preprocess_pdb.slurm, 2_prepare_target_sites.slurm,
# 3_run_masif_mimicry.slurm, and submit_mimicry_pipeline.sh

# --- Host paths ---
MASIF_ROOT="/scratch/ymeng/masif_seed"
MASIF_MIMICRY_ROOT="${MASIF_ROOT}/masif_mimicry"
IMAGE="${MASIF_MIMICRY_ROOT}/masif_mimicry.sif"

# Which pipeline steps to submit (true/false; used by submit_mimicry_pipeline.sh)
RUN_PREPROCESS=false
RUN_TARGET_SITES=false
RUN_MIMICRY_SEARCH=true


# --- Step 1: preprocess target PDB ---
PDB_INPUT="input/021structure.pdb"
PDB_CHAIN="021structure_C_AB"   # p1: chain C; p2: chains A and B
LIGAND="A_021"
LIGAND_SDF="input/A_021.sdf"
TARGET_PREPROCESS_DIR="${MASIF_MIMICRY_ROOT}/data/NUP98"

# --- Step 2: define target sites (paths relative to MASIF_MIMICRY_ROOT for Python CLIs) ---
TARGET_PREPROCESS_DIR_REL="data/NUP98"
TARGET_PDB="021structure_C_AB"
TARGET_CHAIN="C"
TARGET_RESIDUE=728
TARGET_ATOM="CB"
TARGET_PPI_ID="p1"
NUM_POINTS=10
TARGET_SAMPLING_RADIUS=5
TARGET_RUN_DIR_REL="data/NUP98/search_results/021structure_C_AB_"

# --- Step 3: mimicry search ---
DATABASE_DIR="/scratch/ymeng/TED_domainome/output"
DATABASE_SUBSET_DIR="/scratch/ymeng/TED_domainome/output/filtered_intracellular_domainome/subsets"
SEARCH_ARRAY="0-1"

# Desc/dist thresholds (search)
DESC_DIST_CUTOFF=2.0
DESC_DIST_SCORE_CUTOFF=0.45
CA_CLASH_THRESHOLD=1.0
HEAVY_ATOM_CLASH_THRESHOLD=5.0

# --- Derived (do not edit unless needed) ---
MASIF_DB_ROOT="${MASIF_ROOT}/masif"
MASIF_SEED_ROOT="${MASIF_ROOT}/masif_seed_search"
MASIF_DB_SOURCE="${MASIF_DB_ROOT}/source"
MASIF_SEED_SOURCE="${MASIF_SEED_ROOT}/source"
MASIF_MIMICRY_SOURCE="${MASIF_MIMICRY_ROOT}/source"
TARGET_RUN_DIR="${MASIF_MIMICRY_ROOT}/${TARGET_RUN_DIR_REL}"

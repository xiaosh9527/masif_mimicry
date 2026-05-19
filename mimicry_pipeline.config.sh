# Shared configuration for the MaSIF mimicry SLURM pipeline.
# Sourced by 1_preprocess_pdb.slurm, 2_prepare_target_sites.slurm,
# 3_run_masif_mimicry.slurm, 4_postprocess_mimicry.slurm, and submit_mimicry_pipeline.sh

# --- Host paths ---
MASIF_ROOT="/scratch/ymeng/masif_seed"
MASIF_MIMICRY_ROOT="${MASIF_ROOT}/masif_mimicry"
IMAGE="${MASIF_MIMICRY_ROOT}/masif_mimicry.sif"

# Which pipeline steps to submit (true/false; used by submit_mimicry_pipeline.sh)
RUN_PREPROCESS=false
RUN_TARGET_SITES=false
RUN_MIMICRY_SEARCH=true
RUN_POSTPROCESS=true


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
SEARCH_ARRAY="6-10"

# Desc/dist thresholds (search)
DESC_DIST_CUTOFF=2.0
DESC_DIST_SCORE_CUTOFF=0.45
CA_CLASH_THRESHOLD=1.0
HEAVY_ATOM_CLASH_THRESHOLD=5.0

# --- Step 4: post-process search hits (conda env MaSIF; not Apptainer) ---
# POSTPROCESS_TARGET_PDB: chain-C target for clashes/SASA (written by masif_mimicry.define_target_sites)
POSTPROCESS_TARGET_PDB="/scratch/ymeng/masif_seed/masif_mimicry/data/NUP98/search_results/021structure_C_AB_/021structure_AB.pdb"
POSTPROCESS_OUT_DIR="/scratch/ymeng/masif_seed/masif_mimicry/data/NUP98/search_results/postprocess"
POSTPROCESS_CSV_BASE="subset"

# --- Derived (do not edit unless needed) ---
MASIF_DB_ROOT="${MASIF_ROOT}/masif"
MASIF_SEED_ROOT="${MASIF_ROOT}/masif_seed_search"
MASIF_DB_SOURCE="${MASIF_DB_ROOT}/source"
MASIF_SEED_SOURCE="${MASIF_SEED_ROOT}/source"
MASIF_MIMICRY_SOURCE="${MASIF_MIMICRY_ROOT}/source"
TARGET_RUN_DIR="${MASIF_MIMICRY_ROOT}/${TARGET_RUN_DIR_REL}"
POSTPROCESS_OUT_BASENAME="${POSTPROCESS_OUT_DIR}/${POSTPROCESS_CSV_BASE}"

mkdir -p "${POSTPROCESS_OUT_DIR}"
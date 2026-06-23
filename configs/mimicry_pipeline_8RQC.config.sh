# Shared configuration for the MaSIF mimicry SLURM pipelie.
# Sourced by 1_preprocess_pdb.slurm, 2_prepare_target_sites.slurm,
# 3_run_masif_mimicry.slurm, 4_postprocess_mimicry.slurm, and submit_mimicry_pipeline.sh

# --- Host paths ---
MASIF_ROOT="/scratch/ymeng/masif_seed_mimicry"
MASIF_MIMICRY_ROOT="${MASIF_ROOT}/masif_mimicry"
IMAGE="${MASIF_MIMICRY_ROOT}/masif_mimicry.sif"

# Which pipeline steps to submit (true/false; used by submit_mimicry_pipeline.sh)
RUN_PREPROCESS=true
RUN_TARGET_SITES=true
RUN_MIMICRY_SEARCH=true
RUN_POSTPROCESS=true


# --- Step 1: preprocess target PDB ---
PDB_INPUT="data/input/8RQC_AB.pdb"
PDB_CHAIN="8RQC_B_A"   # p1: chain C; p2: chains A and B
LIGAND="A_QFC"
LIGAND_SDF="data/input/8RQC_A_QFC.sdf"
TARGET_PREPROCESS_DIR="${MASIF_MIMICRY_ROOT}/data/8RQC/preprocess"

# --- Step 2: define target sites (paths relative to MASIF_MIMICRY_ROOT for Python CLIs) ---
TARGET_PDB="8RQC_B_A"
TARGET_CHAIN="B"
TARGET_PPI_ID="p1"
NUM_POINTS=15
# Grid selection: residue from any PDB in the same coordinate frame as the preprocessed target
QUERY_PDB="data/input/8RQC_AB.pdb"
QUERY_CHAIN="A"
QUERY_RESIDUE=602
GRID_DISTANCE_CUTOFF=4.0
TARGET_RUN_DIR_REL="data/8RQC/search_results/8RQC_B_A"

# --- Step 3: mimicry search ---
DATABASE_DIR="/work/upthomae/Meng/TED_human_domainome_MaSIF/output"
DATABASE_SUBSET_DIR="/work/upthomae/Meng/TED_human_domainome_MaSIF/output/filtered_intracellular_domainome/subsets"
SEARCH_ARRAY="0-499"

# Desc/dist thresholds (search)
DESC_DIST_CUTOFF=2.0
DESC_DIST_SCORE_CUTOFF=0.45
CA_CLASH_THRESHOLD=1.0
HEAVY_ATOM_CLASH_THRESHOLD=5.0

# --- Step 4: post-process search hits (conda env MaSIF; not Apptainer) ---
# POSTPROCESS_TARGET_PDB: chain-C target for clashes/SASA (written by masif_mimicry.define_target_sites)
# DATABASE_CSV: domain metadata for iface aggregates + merge (filtered to postprocess P1_ids in memory)
DATABASE_CSV="/work/upthomae/Meng/TED_human_domainome_MaSIF/output/filtered_intracellular_domainome/TED_human_domainome_info_intracellular.csv"
POSTPROCESS_TARGET_PDB="${TARGET_RUN_DIR_REL}/8RQC_A.pdb"
POSTPROCESS_OUT_DIR="data/8RQC/postprocess"
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
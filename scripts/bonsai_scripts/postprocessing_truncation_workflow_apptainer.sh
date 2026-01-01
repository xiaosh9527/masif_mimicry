#!/bin/bash

# workflow.sh: Job orchestrator for running MaSIF_seed_mimicry postprocessing 
# and truncation workflow

# 1_postprocessing.sh: launches splits input_pdb_paths.txt into N_ARRAY subsets 
# and submits a job array to calculate postprocessing metrics on matched protein
# models.

# 2_truncate.sh: launches a job array to calculate optimal truncation window using EvoEF2
# and re-calculate metrics on truncated structures.

# -------- Define paths and parameters --------

### Parameters:
LIGAND="B_Y70"   # ligand chain and name
TRUNC_LENGTH=86  # maximum amino acid length of the truncated structures

masif_root=/scratch/shxiao/masif_seed/
IMAGE=$masif_root/masif_mimicry/masif_mimicry.sif

apptainer exec \
  --bind "$masif_root:/workspace" \
  --env LIGAND=$LIGAND \
  --env TRUNC_LENGTH=$TRUNC_LENGTH \
  --pwd /workspace \
  "$IMAGE" /bin/bash -c '
    source /workspace/miniconda/etc/profile.d/conda.sh
    conda activate bonsai

    ### Data paths:
    ROOT=/workspace/masif_mimicry                                             # Root directory
    WORK_DIR=$ROOT/data/template/results/6h0g_C_B_test/P42345-F1-dom-01_A     # IMPORTANT: Working directory containing binder .pdb files to be processed; we use the raw output of poses from P42345 search as an example; needs to be replaced with user-defined directory
    POSTPROCESS_DIR=$WORK_DIR/postprocess                                     # result of postprocessing metrics (default to a subdirectory of the working directory)
    TRUNC_DIR=$WORK_DIR/Truncate_${TRUNC_LENGTH}                              # result of truncation and rescoring metrics (default to a subdirectory of the working directory)
    TARGET_PDB_PATH=$ROOT/data/template/results/6h0g_C_B_test/6h0g_B.pdb      # Target .pdb file MUST INCLUDE THE LIGAND MOLECULE

    ### Script paths: Now defined relative to repo dir
    bonsai_scripts_dir=$ROOT/scripts/bonsai_scripts
    GATHER_SCRIPT=$bonsai_scripts_dir/slurm/gather.sh
    SPLIT_SCRIPT=$bonsai_scripts_dir/slurm/split.sh

    # List the absolute paths of the pdb files in the binder directory as input to postprocessing script (save as input_pdb_paths.txt)
    mkdir -p $POSTPROCESS_DIR
    postprocess_pdb_paths=${POSTPROCESS_DIR}/input_pdb_paths.txt
    find -L "$WORK_DIR" -maxdepth 1 -type f -name "*.pdb" | sort > ${postprocess_pdb_paths}
    echo "PDB paths for postprocessing written to: ${postprocess_pdb_paths}"

    pdb_paths=$(cat ${postprocess_pdb_paths})

    # Run postprocessing and truncation steps
    python $bonsai_scripts_dir/python/postprocess_masif_mimicry.py \
      --pdb_paths $pdb_paths \
      -o $POSTPROCESS_DIR/postprocessed_scores.csv \
      --ligand $LIGAND \
      --target_pdb_path $TARGET_PDB_PATH

    # Truncation and rescoring
    python $bonsai_scripts_dir/python/EvoEF2_truncate.py \
      --clean \
      --input_csv $POSTPROCESS_DIR/postprocessed_scores.csv \
      --length $TRUNC_LENGTH \
      --out_dir $TRUNC_DIR/truncate \
      --temp_dir $TRUNC_DIR/truncate/temp

    # Re-calculate postprocessing metrics on truncated structures
    python $bonsai_scripts_dir/python/proc_trunc_masif_mimicry.py \
      --input $POSTPROCESS_DIR/postprocessed_scores_truncated_results.csv \
      -o $TRUNC_DIR/proc_truncated_scores.csv \
      --ligand $LIGAND
  '
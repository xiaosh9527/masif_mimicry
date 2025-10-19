# MaSIF Mimicry Postprocessing and Truncation Workflow

This repository contains a comprehensive pipeline for postprocessing and analyzing MaSIF surface mimicry search results. 
The workflow processes protein binder structures through multiple stages: postprocessing metrics calculation, optimal truncation, and filtering to identify high-quality binding candidates.

## Overview

The pipeline consists of several interconnected components:

1. **Postprocessing**: Calculate structural and interface metrics for binder structures
2. **Truncation**: Find optimal truncation windows using EvoEF2 folding energy calculations
3. **Filtering**: Apply quality filters based on computed metrics
4. **Sequence retrieval**: Fetch FASTA sequences for final candidates

## Usage Instructions

1. **Setup**: Edit `scripts/config.yaml` with your system paths and executables

2. **Configure Workflow**: Edit parameters in `postprocessing_truncation_workflow.sh`:
   - Set `WORK_DIR`, `BINDER_DIR`, `TARGET_PDB_PATH`
   - Configure `LIGAND` and `TRUNC_LENGTH`
   - Set `DRY_RUN=false` to submit jobs

3. **Run Pipeline**: Submit the workflow script:
   ```bash
   sbatch postprocessing_truncation_workflow.sh
   ```

4. **Filter Results**: Run the R filtering script:
   ```bash
   Rscript R/MaSIF_mimicry_filtering_workflow.R
   ```

5. **Get Sequences**: Fetch FASTA sequences for final candidates:
   ```bash
   python3 fetch_fasta_batch.py \
       --input R/results/proc_trunc_86_6H0F_C_B_250325_filtered.csv \
       --output R/results/proc_trunc_86_6H0F_C_B_250325_filtered_seq.csv
   ```

## Pre-requisites

- **DeepTMHMM**: Pre-computed predictions for human proteome sequences
- **STRIDE**: Executable for secondary structure calculation
- **EvoEF2**: Executable for protein folding free energy calculation
- **SLURM**: Job scheduler for parallel processing
- **R**: For statistical analysis and filtering
- **Python packages**: BioPython, RDKit, pandas, numpy, scikit-image, networkx

## Output Files

- `postprocessed_scores.csv`: Metrics for original binder structures
- `truncated_scores.csv`: Metrics for both original and truncated structures
- `proc_trunc_86_6H0F_C_B_250325_filtered.csv`: Quality-filtered results
- `proc_trunc_86_6H0F_C_B_250325_filtered_seq.csv`: Final results with FASTA sequences 


########################################################
# Details of each component of the workflow:
########################################################

## Workflow Components

### 1. Main Workflow Script: `postprocessing_truncation_workflow.sh`

This is the master script that orchestrates the entire pipeline:

**Key Parameters:**
- `WORK_DIR`: Working directory for storing results
- `BINDER_DIR`: Directory containing MaSIF-matched protein PDB files
- `TARGET_PDB_PATH`: Target protein PDB file
- `LIGAND`: Ligand chain and name (e.g., "B_Y70")
- `TRUNC_LENGTH`: Maximum amino acid length for truncated structures (default: 86)
- `DRY_RUN`: If true, generates SLURM scripts only; if false, submits jobs

**Workflow Steps:**
1. **Postprocessing Phase**: 
   - Generates SLURM job array using `postprocess_wrapper.py`
   - Calculates interface metrics, SASA values, clash counts, and structural properties
   - Outputs: `{WORK_DIR}/postprocess/postprocessed_scores.csv`

2. **Truncation Phase**:
   - Generates SLURM job array using `EvoEF2_trunc_proc_wrapper.py`
   - Finds optimal truncation windows using EvoEF2 folding energy
   - Re-calculates metrics on truncated structures
   - Outputs: `{WORK_DIR}/Truncate_{trunc_length}/truncated_scores.csv`

### 2. Scripts Directory (`/scripts/`)

#### Core Processing Scripts:

**`postprocess_masif_mimicry.py`**: Main postprocessing script
- Calculates clash counts between target and binder structures
- Computes SASA (Solvent Accessible Surface Area) metrics
- Identifies interface residues and calculates interface properties
- Determines intracellular/extracellular localization using DeepTMHMM
- Calculates secondary structure metrics using STRIDE
- Computes geodesic length for structural characterization

**`EvoEF2_truncate.py`**: Truncation optimization script
- Uses EvoEF2 to calculate folding free energies for different truncation windows
- Identifies optimal truncation based on interface residues
- Removes outliers from interface residue lists
- Outputs truncated PDB files with energy scores in filenames

**`proc_trunc_masif_mimicry.py`**: Post-truncation processing script
- Re-calculates all metrics on truncated structures
- Renames metrics with "trunc_" prefix for truncated structures
- Computes difference metrics between original and truncated structures
- Extracts EvoEF2 scores from truncated PDB filenames

#### Wrapper Scripts:

**`postprocess_wrapper.py`**: Generates SLURM job arrays for postprocessing
- Distributes PDB files across array tasks
- Creates parallel processing scripts for scalability

**`EvoEF2_trunc_proc_wrapper.py`**: Generates SLURM job arrays for truncation
- Combines truncation and postprocessing into single job arrays
- Manages intermediate file creation and cleanup

#### Utility Scripts:

**`geodesic_length.py`**: Calculates protein structural elongation
- Uses 3D skeletonization to measure longest dimension
- Computes normalized geodesic length (length/volume^(1/3))

**`ligand_desc_dist_score.py`**: Computes ligand descriptor distance scores
- Measures similarity between query and binder surface descriptors
- Focuses on ligand-binding surface regions

### 3. R Analysis Script: `MaSIF_mimicry_filtering_workflow.R`

Located in `/R/`, this script performs the final filtering step:

**Functions:**
1. **Loads multiple datasets**:
   - Truncated scores (`proc_trunc_86_6H0F_C_B.csv`)
   - Original MaSIF scores (`6h0f_C_B_hp_afdb_merizo_split_search_v2_masif_score.csv`)
   - Filtering thresholds (`trunc_filters_mimicry_6H0F_C_B_250325.csv`)

2. **Calculates additional metrics**:
   - Interface-to-terminus distances (minimum and mean)
   - UniProt exclusion flags
   - Domain exclusion flags

3. **Applies quality filters**:
   - Filters based on manually determined thresholds
   - Combines MaSIF scores with postprocessing metrics
   - Outputs filtered results

### 4. Sequence Retrieval: `fetch_fasta_batch.py`

**Purpose**: Retrieves FASTA sequences for filtered protein candidates

**Features:**
- Batch processing of UniProt IDs using UniProt REST API
- Handles large datasets efficiently with batching (50 IDs per request)
- Extracts truncated sequences based on truncation coordinates
- Robust error handling for API failures

**Usage:**
```bash
python3 fetch_fasta_batch.py \
    --input /path/to/filtered_results.csv \
    --output /path/to/results_with_sequences.csv
```

## Configuration

### `config.yaml`
Contains system-specific paths and executable locations:
- `deeptmhmm_dir`: Path to pre-computed DeepTMHMM predictions
- `python_path`: Python executable path
- `stride_exec`: STRIDE executable command
- `EvoEF2_dir`: EvoEF2 installation directory
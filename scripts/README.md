1. Postprocess_masif_mimicry.py:
• Postprocessing scores were calculated for the resulting matched proteins
• query_preprocessed_root in this script is a vestige from adopting from MaSIF neosurf workflow. It is a mandatory argument but is not used for actual calculation of metrics.
• deepTMHMM reference table generation script: see hp_list_DeepTMHMM folder:
• prep_fasta: write {uniprot_id}.fasta for all accessions in hp_list.txt to /fasta folder
• run_DeepTMHMM_array.sh requires fasta dir, subset dir (for job array), and deeptmhmm package
2. 4_EvoEF2_truncate.py:
• Truncate by sliding window, and calculate folding free energy of each truncation by EvoEF2.
• 4_fix_broken_auth_chain.py: The .pdb chains are broken in the domainome database, so inside 4_EvoEF2_truncate.py I also call subprocess run to fix the broken chain ids first, before truncating using the fixed .pdb as input for EvoEF2
• .sh bash script to read postprocessed score.csv, split into subsets, run truncation, read truncation output, and write to subset_*.csv
3. Proc_trunc_masif_mimicry.py:
• Re-calculate some of the postprocessing metrics based on the truncated structures.
4. Filters:
• The filtering metrics cutoffs are determined manually by human inspection of the structures to decide what the cutoffs should be to filter out binder models that makes no structural sense. The final range/thresholds used are recorded in .csv files. Note the .csv file contains the range of all the metrics that had been calculated, even if they were not used to filter out any hit (then the min and max would be set to the maximum range).
• After filtering, the 2000 truncated binders with the highest MaSIF.score were selected for sequence extraction.
5. Fetch_fasta_batch.py:
• Use uniprot_id, trunc_resi_start, and trunc_resi_end columns in the .csv file written by proc_trunc_masif_mimicry.py to construct batch PDB API call to retrieve their protein sequences. Depreciated uniprot_id will return NA, in which case I manually found the protein sequences.

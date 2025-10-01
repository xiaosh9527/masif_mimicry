import os, subprocess, requests, argparse

import pandas as pd
import numpy as np

from Bio.PDB import *

def create_parser():
    parser = argparse.ArgumentParser(description='Split AlphaFold models into domains based on PAE.')
    parser.add_argument('--uniprot_id', type=str, required=True, help='UniProt ID of the protein to process.')
    parser.add_argument('--pae_cutoff', type=float, default=15.0, help='PAE cutoff to define domains.')
    parser.add_argument('--input_dir', type=str, default='./input/', help='Directory to store input files.')
    return parser

def retrive_data_from_url(url, outfile):
    response = requests.get(url)
    if response.status_code == 200:
        with open(outfile, 'wb') as f:
            f.write(response.content)
    else:
        raise Exception(f"Failed to retrieve data from {url}. Status code: {response.status_code}")

def download_af_db_data(uniprot_id, out_pdb_folder='./', out_pae_folder='./'):
    af_pdb_query = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb'
    af_pae_query = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-predicted_aligned_error_v4.json'


    out_pdb_name = os.path.join(out_pdb_folder, f"AF-{uniprot_id}-F1-model_v4.pdb")
    out_pae_name = os.path.join(out_pae_folder, f"AF-{uniprot_id}-F1-predicted_aligned_error_v4.json")

    retrive_data_from_url(af_pdb_query, out_pdb_name)
    retrive_data_from_url(af_pae_query, out_pae_name)
    
def row_to_list(row):
    # Converts row df to list
    return [int(r) for r in row.dropna().tolist()]

def expand_cluster(selected_resnums, max_neighbor_dist=5):
    expanded_selection = []
    for i,res in enumerate(selected_resnums):
        expanded_selection.append(res)
        try:
            next_res = selected_resnums[i + 1]
            neighbor_dist = next_res - res
            if neighbor_dist == 1:
                continue # Consecutive resnums
            elif neighbor_dist <= max_neighbor_dist:
                missing_neighbors = [n for n in range(res + 1, next_res)]
                expanded_selection.extend(missing_neighbors)
        except IndexError:
            break
    return expanded_selection

def extract_residues(input_pdb, residue_list, output_pdb, cap_terminis=False):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("input_structure", input_pdb)
    
    # To get the last residue number
    last_chain = structure[0].child_list[-1]
    last_residue_number = last_chain.child_list[-1].id[1]
    
    # Create a new Structure object to store the selected residues
    output_structure = Structure.Structure("output_structure")
    
    if residue_list == 'all':
        residue_list = [residue.id[1] for residue in structure.get_residues()]
        # print(residue_list)
        
    # Iterate over all models in the input structure
    for model in structure:
        # Create a new Model object for the output structure
        output_model = Model.Model(model.id)
        
        # Iterate over all chains in the model
        for chain in model:
            # Create a new Chain object for the output structure
            output_chain = Chain.Chain(chain.id)
                        
            # Iterate over all residues in the chain
            all_residues_in_full_pdb = [residue for m in structure for chain in model for residue in chain]
            for residue_index, residue in enumerate(all_residues_in_full_pdb):
                # Check if the residue number is in the list
                if residue.id[1] in residue_list:
                    inx_in_residue_list = residue_list.index(residue.id[1])
                    
                    if cap_terminis:
                        # Check if this is an N-terminal residue
                        if residue.id[1] == 1: # Real N-terminal
                            pass
                        elif inx_in_residue_list == 0: # Non-real N-terminal
                            #print(f'First fragment residue. Non-real N terminal frag {residue}')
                            new_ter_residue = Residue.Residue((' ', residue.id[1] - 1, ' '), 'ACE', residue.id[1] -1)
                            for atom in all_residues_in_full_pdb[residue_index - 1]:
                                if atom.id in ['C', 'N', 'O', 'CA']:
                                    if atom.id == 'CA':
                                        atom.id = 'CH3' # Change name to coincide with ACE atomnames
                                    # Assign a unique serial number to each atom
                                    atom.set_serial_number(len(new_ter_residue) + 1)
                                    new_ter_residue.add(atom.copy())
                            # Add the new residue to the output chain
                            output_chain.add(new_ter_residue)
                        elif residue.id[1] -1 != residue_list[inx_in_residue_list - 1]: # Non-real N-terminal
                            #print(f'Non-real N terminal frag {residue}')
                            new_ter_residue = Residue.Residue((' ', residue.id[1] - 1, ' '), 'ACE', residue.id[1] -1)
                            for atom in all_residues_in_full_pdb[residue_index - 1]:
                                if atom.id in ['C', 'N', 'O', 'CA']:
                                    if atom.id == 'CA':
                                        atom.id = 'CH3' # Change name to coincide with ACE atomnames
                                    # Assign a unique serial number to each atom
                                    atom.set_serial_number(len(new_ter_residue) + 1)
                                    new_ter_residue.add(atom.copy())
                            # Add the new residue to the output chain
                            output_chain.add(new_ter_residue)
                            
                    # Create a new residue with unique serial numbers for atoms
                    new_residue = Residue.Residue((' ', residue.id[1], ' '), residue.get_resname(), residue.id[1])
                    for atom in residue:
                        # Assign a unique serial number to each atom
                        atom.set_serial_number(len(new_residue) + 1)
                        new_residue.add(atom.copy())
                    # Add the new residue to the output chain
                    output_chain.add(new_residue)
                    
                    if cap_terminis:
                        # Check if this is a terminal residue
                        if residue.id[1] == last_residue_number: # Real C-terminal
                            pass
                        elif inx_in_residue_list == len(residue_list) -1: # Non-real C-terminal
                            new_ter_residue = Residue.Residue((' ', residue.id[1] + 1, ' '), 'NME', residue.id[1] + 1)
                            for atom in all_residues_in_full_pdb[residue_index + 1]:
                                if atom.id in ['N', 'CA']:
                                    if atom.id == 'CA':
                                        atom.id = 'CH3' # Change name to coincide with ACE atomnames
                                    # Assign a unique serial number to each atom
                                    atom.set_serial_number(len(new_ter_residue) + 1)
                                    new_ter_residue.add(atom.copy())
                            # Add the new residue to the output chain
                            output_chain.add(new_ter_residue)

                        elif residue.id[1] + 1 != residue_list[inx_in_residue_list + 1]: # Non-real C-terminal
                            #print(f'Non-real C terminal frag {residue}')
                            new_ter_residue = Residue.Residue((' ', residue.id[1] + 1, ' '), 'NME', residue.id[1] + 1)
                            for atom in all_residues_in_full_pdb[residue_index + 1]:
                                if atom.id in ['N', 'CA']:
                                    if atom.id == 'CA':
                                        atom.id = 'CH3' # Change name to coincide with ACE atomnames
                                    # Assign a unique serial number to each atom
                                    atom.set_serial_number(len(new_ter_residue) + 1)
                                    new_ter_residue.add(atom.copy())
                            # Add the new residue to the output chain
                            output_chain.add(new_ter_residue)
                    
            
            # Add the output chain to the output model
            output_model.add(output_chain)
        
        # Add the output model to the output structure
        output_structure.add(output_model)

    # Write out the selected residues to a new PDB file
    io = PDBIO()
    io.set_structure(output_structure)
    io.save(output_pdb)
    
def truncate_low_plddt_loops(input_path, output_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('', input_path)
    dssp = DSSP(structure[0], input_path, dssp="/work/lpdi/bin/sequence/dssp")
    
    dssp_label = [x[2] for x in dssp]
    plddt_list = [atom.get_bfactor() for atom in structure.get_atoms() if atom.get_name() == 'CA']
    all_atoms = [atom for atom in structure.get_atoms()]
    
    # First filter: remove low confidence regions from N- and C- terminus
    del_res_list = np.zeros((len(plddt_list)), dtype=bool)
    for i, (sec_label, plddt) in enumerate(zip(dssp_label, plddt_list)):
        if plddt > 70.0 or sec_label not in ['-', 'G', 'S', 'T']:
            break
        else:
            del_res_list[i] = True
            
    for i, (sec_label, plddt) in reversed(list(enumerate(zip(dssp_label, plddt_list)))):
        if plddt > 70.0 or sec_label not in ['-', 'G', 'S', 'T']:
            break
        else:
            del_res_list[i] = True
    
    all_res = np.array([res.get_id()[1] for res in structure.get_residues()])
    remain_res = all_res[np.where(del_res_list==False)[0]]
    
    # Second filter: if remaining residue length is shorter than 5, remove from the list
    if len(remain_res) >= 5:
        # remain_avg_plddt = np.array(plddt_list)[np.where(del_res_list==False)[0]].mean()
        io = PDBIO()
        for atom in Selection.unfold_entities(structure, 'A'):
            if atom.get_parent().get_id()[1] not in remain_res:
                parent = atom.get_parent()
                parent.detach_child(atom.get_id())
        io.set_structure(structure)
        io.save(output_path)
        # print(f'{pdb_file} removed due to low avg plddt!')
        return True
    else:
        # print(f'{pdb_file} removed!')
        return False

def main(uniprot_id: str = 'P0DTC2', pae_cutoff: float = 15.0):
    download_af_db_data(uniprot_id, out_pae_folder=PAE_FOLDER, out_pdb_folder=RAW_PDB_FOLDER)

    pdb_file = os.path.join(RAW_PDB_FOLDER, f'AF-{uniprot_id}-F1-model_v4.pdb')
    pae_file = os.path.join(PAE_FOLDER, f'AF-{uniprot_id}-F1-predicted_aligned_error_v4.json')
    tmp_file = os.path.join(TMP_FOLDER, f'AF-{uniprot_id}-F1-predicted_aligned_clusters.csv')

    subprocess.run(['python', '/work/lpdi/users/shxiao/pae_to_domains/pae_to_domains.py', pae_file, '--pae_cutoff', str(pae_cutoff), '--output_file', tmp_file])

    df_clust = pd.read_csv(tmp_file, header=None)
    list_of_clust = df_clust.apply(lambda row: row_to_list(row), axis=1)

    for i, clust in enumerate(list_of_clust):
        print(f'Step 1: Processing cluster {i+1} for {uniprot_id}...')
        out_file_1 = os.path.join(TMP_FOLDER, f'{uniprot_id}-F1-dom{i+1}_tmp1.pdb')
        expanded_clust = expand_cluster(clust, max_neighbor_dist=5)
        extract_residues(pdb_file, expanded_clust, out_file_1, cap_terminis=False)
        print(f'Step 2: Truncating low pLDDT loops for cluster {i+1}...')
        out_file_2 = os.path.join(TMP_FOLDER, f'{uniprot_id}-F1-dom{i+1}_tmp2.pdb')
        done = truncate_low_plddt_loops(out_file_1, out_file_2)
        if done:
            print(f'Step 3: Adding missing terminal atoms NME & ACE to N & C terminals...')
            out_file_3 = os.path.join(FRAG_FOLDER, f'{uniprot_id}-F1-dom{i+1}.pdb')
            extract_residues(out_file_2, 'all', out_file_3, cap_terminis=True)
        else:
            print(f'Step 3 skipped due for disorganized regions!')
            pass

if __name__ == '__main__':        
    parser = create_parser()
    args = parser.parse_args()  

    TMP_FOLDER = os.getenv('TMPDIR')

    RAW_PDB_FOLDER = f'{args.input_dir}/raw_pdb/'
    PAE_FOLDER = f'{args.input_dir}/pae_data/'    
    FRAG_FOLDER = f'{args.input_dir}/fragments/'

    os.makedirs(RAW_PDB_FOLDER, exist_ok=True)
    os.makedirs(PAE_FOLDER, exist_ok=True)
    os.makedirs(FRAG_FOLDER, exist_ok=True)

    main(uniprot_id=args.uniprot_id, pae_cutoff=args.pae_cutoff)
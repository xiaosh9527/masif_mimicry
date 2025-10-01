import pandas as pd
import requests
import time

def fetch_batch_fasta(uniprot_ids):
    """
    Fetches FASTA sequences for a list of UniProt IDs in batch mode using the UniProt REST API.
    The query is built as accession:(ID1 OR ID2 OR ...).
    """
    # Build the query string: e.g. accession:(P12345 OR Q8N158 OR O14920)
    query_ids = " OR ".join(uniprot_ids)
    query = f"accession:({query_ids})"
    url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query={query}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        print(f"Error fetching batch for IDs {uniprot_ids}: HTTP {response.status_code}")
        return None

def parse_fasta(fasta_str):
    """
    Parses a FASTA-formatted string and returns a dictionary mapping UniProt ID to sequence.
    Assumes that the header is formatted as >sp|P12345|... and extracts the second field as the ID.
    """
    sequences = {}
    current_id = None
    seq_lines = []
    for line in fasta_str.splitlines():
        if line.startswith(">"):
            if current_id:
                sequences[current_id] = "".join(seq_lines)
            # Extract UniProt ID (assumes header like >sp|P12345|...)
            parts = line.split("|")
            if len(parts) >= 2:
                current_id = parts[1]
            else:
                # Fallback: take the first word after '>'
                current_id = line[1:].split()[0]
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if current_id:
        sequences[current_id] = "".join(seq_lines)
    return sequences

def main():
    # Read the input CSV file (adjust filename as needed)
    df = pd.read_csv("/Users/ymeng/Desktop/ThomaLab/R_projects/MaSIF_surface_mimicry_inspection/results/proc_trunc_86_6H0F_C_B_250325_pending.csv")
    
    # Get unique UniProt IDs to minimize duplicate requests
    unique_ids = list(df["uniprot_id"].unique())
    
    # We'll break the list into smaller batches to avoid excessively long URLs
    batch_size = 50
    sequences = {}
    for i in range(0, len(unique_ids), batch_size):
        batch_ids = unique_ids[i:i+batch_size]
        fasta_data = fetch_batch_fasta(batch_ids)
        if fasta_data:
            batch_sequences = parse_fasta(fasta_data)
            sequences.update(batch_sequences)
        else:
            print(f"Batch starting at index {i} failed to retrieve data.")
        time.sleep(1)  # Sleep 1 second between batches to be polite to the server
    
    # Create a new column for the truncated sequence
    truncated_seqs = []
    for index, row in df.iterrows():
        uniprot_id = row["uniprot_id"]
        start = int(row["trunc_resi_start"])
        end = int(row["trunc_resi_end"])
        full_seq = sequences.get(uniprot_id)
        if full_seq:
            # Adjust for 0-based indexing in Python: subtract 1 from start
            trunc_seq = full_seq[start - 1:end]
        else:
            trunc_seq = None
        truncated_seqs.append(trunc_seq)
    
    df["turnc_sequence"] = truncated_seqs
    df.to_csv("/Users/ymeng/Desktop/ThomaLab/R_projects/MaSIF_surface_mimicry_inspection/results/proc_trunc_86_6H0F_C_B_250325_pending_seq.csv", index=False)
    print("Updated CSV file saved as 'output.csv'.")

if __name__ == "__main__":
    main()
import subprocess
from typing import List



def get_stride_annotations(
        pdb_path: str,
        stride_exec: str = 'stride',
    ):
    """
    Annotations:
    H 	    Alpha helix
    G 	    3-10 helix
    I 	    PI-helix
    E 	    Extended conformation (beta strand)
    B or b 	Isolated bridge
    T 	    Turn
    C 	    Coil (none of the above)
    """

    result = subprocess.run(
        [stride_exec, pdb_path], 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE, 
        encoding='utf-8', 
        check=True,
    )

    # Build a mapping from residue number to SS assignment from lines starting with "ASG"
    stride_mapping = {}
    for line in result.stdout.splitlines():
        if line.startswith("ASG"):
            parts = line.split()
            res_num = int(parts[3])
            ss_code = parts[5]  # parts[5] is the secondary structure code from stride output.
            assert len(ss_code) == 1
            stride_mapping[res_num] = ss_code

    return stride_mapping


def get_mdtraj_annotation(pdb_path: str) -> list:
    import mdtraj as md
    """
    Annotations:
    H: Helix (H, G or I in the original DSSP code)
    E: Strand (E or B in the original DSSP code)
    C: Coil (T or S in the original DSSP code)
    NA: Assigned if it's not a protein residue
    """
    traj = md.load(pdb_path)  # use first model; use md.load_pdb for speed if no altlocs
    dssp_codes = md.compute_dssp(traj, simplified=True)  # array shape (n_frames, n_res)
    dssp_codes = dssp_codes[0]  # first frame
    return dssp_codes


def find_sse(resi: List[int], annotation: List[str]) -> List[dict]:
    """Find contiguous segments with the same secondary structure."""
    
    segments = []
    for i, label in zip(resi, annotation):
        # previous segment continues
        if (len(segments) > 0) and (i == segments[-1]['end'] + 1 and label == segments[-1]['label']):
            segments[-1]['end'] = i

        # new segment 
        else:
            segments.append({'start': i, 'end': i, 'label': label})

    return segments

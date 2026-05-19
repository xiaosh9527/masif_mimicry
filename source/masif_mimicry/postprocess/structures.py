"""Shared structure I/O for mimicry post-processing."""

from pathlib import Path
from typing import Union

from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.StructureBuilder import StructureBuilder

STANDARD_AA = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
})


def filter_structure_to_standard_aa(struct, standard_aa=None):
    """
    Return a new Structure with only standard amino acid residues.
    Preserves chain IDs and residue numbering.
    """
    if standard_aa is None:
        standard_aa = STANDARD_AA

    model = struct[0] if isinstance(struct, Structure) else struct
    builder = StructureBuilder()
    builder.init_structure("filtered")
    builder.init_model(0)

    for chain in model.get_chains():
        builder.init_chain(chain.id)

    out_struct = builder.get_structure()
    out_model = out_struct[0]

    for chain in model.get_chains():
        chain_id = chain.id
        for res in chain.get_residues():
            if res.get_resname() in standard_aa:
                out_model[chain_id].add(res.copy())

    return out_struct


def maybe_load_structure(path_or_structure: Union[str, Path, Structure]):
    if isinstance(path_or_structure, Structure):
        return path_or_structure
    return PDBParser(QUIET=True).get_structure("", str(path_or_structure))

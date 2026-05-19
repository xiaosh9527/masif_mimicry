import os

from Bio.PDB import PDBParser


def parse_partner_chain_ids(pdb_identifier, target_ppi_id):
    """
    Return partner chain id(s) for clash counting from a MaSIF PDB identifier.
    """
    parts = pdb_identifier.split("_")
    if target_ppi_id == "p1":
        if len(parts) < 3:
            raise ValueError(
                f"Cannot parse partner chains from {pdb_identifier!r} with target_ppi_id p1"
            )
        raw = parts[2]
    else:
        if len(parts) < 2:
            raise ValueError(
                f"Cannot parse partner chains from {pdb_identifier!r} with target_ppi_id p2"
            )
        raw = parts[1]
    return _chain_ids_for_filter(raw)


def _chain_ids_for_filter(target_chain):
    """Normalize chain argument(s) to a list of chain ids for membership tests."""
    if isinstance(target_chain, (list, tuple)):
        return [str(c) for c in target_chain]
    if isinstance(target_chain, str):
        if len(target_chain) == 1:
            return [target_chain]
        return list(target_chain)
    return [str(target_chain)]


def get_filtered_target_structure(target_pdb_path, target_chain, cache):
    """Parse target PDB once and cache chain-filtered structure for clash counting."""
    chain_ids = _chain_ids_for_filter(target_chain)
    key = (os.path.abspath(target_pdb_path), tuple(chain_ids))
    if key not in cache:
        pdb_parser = PDBParser(QUIET=True)
        target_structure = pdb_parser.get_structure("", target_pdb_path)
        chains_to_remove = [
            chain for chain in target_structure.get_chains()
            if chain.id not in chain_ids
        ]
        for chain in chains_to_remove:
            target_structure[0].detach_child(chain.id)
        cache[key] = target_structure
    return cache[key]

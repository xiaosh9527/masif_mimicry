from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.StructureBuilder import StructureBuilder

from masif_mimicry.postprocess.structures import maybe_load_structure


def get_sasa(structure, sub_surface_key=None):
    for res in structure.get_atoms():
        computed = hasattr(res, 'sasa')
        break
    if not computed:
        ShrakeRupley().compute(structure, level="A")
    if sub_surface_key is not None:
        structure = structure[sub_surface_key]
    sasa = sum(a.sasa for a in structure.get_atoms())
    return sasa


def extract_ligand_residue(pdb_struct, ligand_chain, ligand_name):
    ligand = [res for res in pdb_struct[ligand_chain].get_residues() if res.get_resname() == ligand_name]
    assert len(ligand) == 1
    return ligand[0]


def extract_ligand(pdb_struct, ligand_chain, ligand_name):
    builder = StructureBuilder()
    builder.init_structure('ligand')
    builder.init_model(0)

    builder.init_chain("L")  # for ligand

    out_struct = builder.get_structure()[0]
    ligand = [res for res in pdb_struct[ligand_chain].get_residues() if res.get_resname() == ligand_name]
    assert len(ligand) == 1
    copy_residues(ligand, out_struct, output_chain="L")

    return out_struct


def copy_residues(input_residues, output_struct, output_chain):

    if not isinstance(input_residues, list):
        input_residues = input_residues.get_residues()

    for i, res in enumerate(input_residues):
        _res = res.copy()
        _res.id = (" ", i + 1, " ")
        output_struct[output_chain].add(_res)


def merge_structures(struct_1, struct_2=None, key1="T", key2="B", ligand=None, ligand_def=None):
    builder = StructureBuilder()
    builder.init_structure('complex')
    builder.init_model(0)

    builder.init_chain(key1)  # for target
    if struct_2 is not None:
        builder.init_chain(key2)  # for binder
    if ligand is not None or ligand_def is not None:
        builder.init_chain("L")  # for ligand

    out_struct = builder.get_structure()[0]

    if ligand_def is not None:
        # remove ligand from target structure to not count it twice
        struct_1 = [res for res in struct_1.get_residues() if not (res.get_resname() == ligand_def["name"] and res.parent.id == ligand_def["chain"])]

    copy_residues(struct_1, out_struct, output_chain=key1)
    if struct_2 is not None:
        copy_residues(struct_2, out_struct, output_chain=key2)
    if ligand is not None:
        copy_residues(ligand, out_struct, output_chain="L")

    return out_struct


def get_isolated_chain(target_struct):
    builder = StructureBuilder()
    builder.init_structure('target')
    builder.init_model(0)

    builder.init_chain("T")  # for target
    out_struct = builder.get_structure()[0]

    copy_residues(target_struct, out_struct, output_chain="T")

    return out_struct


def compute_sasa_values(target_path, matched_path, ligand_def):
    target_ligand = maybe_load_structure(target_path)[0]
    ligand = extract_ligand(target_ligand, ligand_def["chain"], ligand_def["name"])
    binder = maybe_load_structure(matched_path)[0]
    whole_complex = merge_structures(struct_1=target_ligand, struct_2=binder, ligand=ligand, ligand_def=ligand_def)
    target = get_isolated_chain(whole_complex["T"])
    binder_ligand = merge_structures(struct_1=binder, struct_2=None, ligand=ligand, key1="B")
    binder_target = merge_structures(struct_1=target_ligand, struct_2=binder, key1="T", ligand_def=ligand_def)
    target_ligand = merge_structures(struct_1=target_ligand, struct_2=None, ligand=ligand, ligand_def=ligand_def)

    # all possible full surface areas
    t_l_sasa = get_sasa(target_ligand)
    b_l_sasa = get_sasa(binder_ligand)
    t_sasa = get_sasa(target)
    b_sasa = get_sasa(binder)
    l_sasa = get_sasa(ligand)
    t_b_sasa = get_sasa(binder_target)
    complex_sasa = get_sasa(whole_complex)

    # surface areas contributing to different complexes
    t_in_complex_sasa = get_sasa(whole_complex, sub_surface_key="T")
    b_in_complex_sasa = get_sasa(whole_complex, sub_surface_key="B")
    l_in_complex_sasa = get_sasa(whole_complex, sub_surface_key="L")
    t_in_tb_sasa = get_sasa(binder_target, sub_surface_key="T")
    b_in_tb_sasa = get_sasa(binder_target, sub_surface_key="B")
    l_in_lb_sasa = get_sasa(binder_ligand, sub_surface_key="L")
    b_in_lb_sasa = get_sasa(binder_ligand, sub_surface_key="B")
    l_in_tl_sasa = get_sasa(target_ligand, sub_surface_key="L")
    t_in_tl_sasa = get_sasa(target_ligand, sub_surface_key="T")

    # buried single surface areas in the complex
    buried_t = t_sasa - t_in_complex_sasa
    buried_b = b_sasa - b_in_complex_sasa
    buried_l = l_sasa - l_in_complex_sasa
    buried_total = buried_t + buried_b + buried_l

    # buried single surface area in binary
    buried_t_tb = t_sasa - t_in_tb_sasa
    buried_b_tb = b_sasa - b_in_tb_sasa
    buried_l_lb = l_sasa - l_in_lb_sasa
    buried_b_lb = b_sasa - b_in_lb_sasa
    buried_l_tl = l_sasa - l_in_tl_sasa
    buried_t_tl = t_sasa - t_in_tl_sasa

    # lig buried
    buried_l_in_tl_by_complex = l_in_tl_sasa - l_in_complex_sasa
    buried_tl_by_complex = t_l_sasa - t_in_complex_sasa

    # ligand contribution to different interface areas
    ligand_contrib_total = buried_l / buried_total
    lig_contrib_iface_to_binder = buried_l_in_tl_by_complex / buried_tl_by_complex

    delta_sasa = (t_l_sasa + b_l_sasa) - complex_sasa

    out = {
        # ---surface areas unary---
        "target_unbound_sasa": t_sasa,
        "binder_unbound_sasa": b_sasa,
        "ligand_unbound_sasa": l_sasa,
        # ---surface areas binary---
        "target_ligand_sasa": t_l_sasa,
        "binder_ligand_sasa": b_l_sasa,
        # "target_binder_sasa": t_b_sasa,
        # ---surface area ternary---
        "complex_sasa": complex_sasa,
        # ---sub surface areas of ternary complex---
        "target_in_complex_sasa": t_in_complex_sasa,
        "binder_in_complex_sasa": b_in_complex_sasa,
        "ligand_in_complex_sasa": l_in_complex_sasa,
        # ---sub surface areas of binary complexes---
        # "target_in_tb_sasa": t_in_tb_sasa,
        # "binder_in_tb_sasa": b_in_tb_sasa,
        "ligand_in_lb_sasa": l_in_lb_sasa,
        "binder_in_lb_sasa": b_in_lb_sasa,
        "ligand_in_tl_sasa": l_in_tl_sasa,
        "target_in_tl_sasa": t_in_tl_sasa,
        # ---buried surface areas---
        "target_buried_in_complex": buried_t,
        "binder_buried_in_complex": buried_b,
        "ligand_buried_in_complex": buried_l,
        # "buried_total": buried_total,
        "target_buried_in_tb": buried_t_tb,
        "binder_buried_in_tb": buried_b_tb,
        "ligand_buried_in_lb": buried_l_lb,
        "binder_buried_in_lb": buried_b_lb,
        "ligand_buried_in_tl": buried_l_tl,
        "target_buried_in_tl": buried_t_tl,
        # ---ligand contributions to interfaces---
        # "ligand_iface_contribution_total": ligand_contrib_total,
        "ligand_iface_contribution": lig_contrib_iface_to_binder,
        # ---dSASA---
        "delta_sasa": delta_sasa
    }
    return out

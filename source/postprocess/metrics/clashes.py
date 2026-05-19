import numpy as np
from rdkit import Chem


def coord_and_radii(rdmol, ignore={'H'}):
    _periodic_table = Chem.GetPeriodicTable()

    coord = rdmol.GetConformer().GetPositions()
    radii = np.array([_periodic_table.GetRvdw(a.GetSymbol()) for a in rdmol.GetAtoms()])

    mask = np.array([a.GetSymbol() not in ignore for a in rdmol.GetAtoms()])
    coord = coord[mask]
    radii = radii[mask]

    assert coord.shape[0] == radii.shape[0]
    return coord, radii


def count_clashes(pdb1, pdb2, strictness=1.0):
    """Compute number of clashing heavy atom pairs."""

    protein1 = Chem.MolFromPDBFile(str(pdb1), sanitize=False)
    protein2 = Chem.MolFromPDBFile(str(pdb2), sanitize=False)

    coord1, radii1 = coord_and_radii(protein1)
    coord2, radii2 = coord_and_radii(protein2)

    dist = np.sqrt(np.sum((coord1[:, None, :] - coord2[None, :, :]) ** 2, axis=-1))

    clashes = dist < strictness * (radii1[:, None] + radii2[None, :])
    num_clashes = np.sum(clashes)

    return num_clashes
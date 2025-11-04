import tempfile

from biotite.structure.io.pdb import PDBFile
from biotite.structure import BondList
from rdkit import Chem

import hydride
import numpy as np


def prepare_ligand_from_protonated_rdkit_mol(
    rdkit_mol,
    out_pdb="./tmp/ligand_pH7_h_hydride.pdb"
):
    """
    Function for adding hydrogens to ligand molecule given by SMILES format.
    Achieved with hydride.
    """

    print("---- Adding hydrogens with hydride ----")

    # saving ligand without Hs as PDB
    tmp_pdb = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    print(tmp_pdb)
    Chem.MolToPDBFile(rdkit_mol, tmp_pdb)

    # loading PDB in Biotite
    pdb = PDBFile.read(tmp_pdb)
    atom_array = pdb.get_structure(model=1)

    # creating BondList for hydride
    n_atoms = atom_array.array_length()
    bond_list = BondList(n_atoms)

    for bond in rdkit_mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        order = bond.GetBondTypeAsDouble()
        bond_list.add_bond(a1, a2, order)

    atom_array.bonds = bond_list

    # loading formal charges from RDKit - charges still appears to be zeros
    rdkit_atoms = list(rdkit_mol.GetAtoms())
    charges = np.zeros(len(atom_array), dtype=float)
    min_len = min(len(rdkit_atoms), len(atom_array))
    for i in range(min_len):
        charges[i] = float(rdkit_atoms[i].GetFormalCharge())
    atom_array.charge = charges

    # adding Hs with hydride
    try:
        print("---- Starting hydride ----")
        mask = [True] * len(atom_array)  # this ensures all atoms gets Hs?
        hydrogenated_array, mask_out = hydride.add_hydrogen(
            atom_array, mask=mask)

        # copying formal charges to new molecule with Hs
        new_charges = np.zeros(len(hydrogenated_array), dtype=float)
        new_charges[:len(atom_array)] = atom_array.charge
        hydrogenated_array.charge = new_charges

        # saving result PDB
        pdb.set_structure(hydrogenated_array)
        pdb.write(out_pdb)
        print(f"---- Ligand with Hs saved to {out_pdb} ----")

        return out_pdb

    except Exception as e:
        print(f"Hydride failed: {e}")
        raise

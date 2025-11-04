from rdkit import Chem
from rdkit.Chem import AllChem

from dimorphite_dl import protonate_smiles


def protonate_with_dimorphite(path: str, pH_min: float = 7.0, pH_max: float = 7.0, max_variants: int = 5):
    """
    Function for loading and protonating molecule given by path to PDB file.
    Achieved with dimorphite_dl.
    """
    # loading ligand from PDB file (without H)
    orig = Chem.MolFromPDBFile(path, removeHs=True, sanitize=False)
    if orig is None:
        raise ValueError(f"RDKit failed to load PDB file: {path}")

    # removing H and generating SMILES
    orig_no_h = Chem.RemoveAllHs(orig)
    base_smiles = Chem.MolToSmiles(orig_no_h, isomericSmiles=True)
    print(f"---- Base SMILES: {base_smiles} ----")

    # protonating SMILES with dimorphite
    prot_list = protonate_smiles(base_smiles,
                                 ph_min=pH_min,
                                 ph_max=pH_max,
                                 label_states=True,
                                 max_variants=max_variants)
    if not prot_list:
        raise RuntimeError(
            "Dimorphite-DL failed to protonate SMILES")

    # protonated SMILES (choosing from variants, taking first variant)
    protonated_smiles = prot_list[0]
    print(
        f"---- Protonated SMILES (pH {pH_min, pH_max}): {protonated_smiles} ----")

    # creating Mol from SMILES
    prot_mol_no_h = Chem.MolFromSmiles(protonated_smiles)
    if prot_mol_no_h is None:
        raise RuntimeError(
            "Failed to convert from SMILES to Mol")

    # sanitizing
    try:
        Chem.SanitizeMol(prot_mol_no_h)
    except Exception:
        print("Sanitization failed")

    # removing Hs
    # prot_mol_no_h = Chem.RemoveAllHs(prot_mol)
    # print(prot_mol_no_h)

    # 3D conformer - without this ligand has 9 Hs, with this ligand has 10 Hs
    if prot_mol_no_h.GetNumConformers() == 0:
        prot_mol_no_h = Chem.AddHs(prot_mol_no_h, addCoords=False)
        AllChem.EmbedMolecule(prot_mol_no_h, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(prot_mol_no_h)
        prot_mol_no_h = Chem.RemoveAllHs(prot_mol_no_h)

    return prot_mol_no_h

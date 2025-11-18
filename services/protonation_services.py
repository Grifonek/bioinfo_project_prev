from rdkit import Chem
from dimorphite_dl import protonate_smiles

class ProtonationService:
    """
    Service that protonates ligands using dimorphite_dl.
    Methods return an RDKit Mol without explicit hydrogens,
    leaving the addition of hydrogens to the HydrideService.
    """

    def __init__(self, max_variants: int = 5):
        self.max_variants = max_variants

    def protonate_with_dimorphite(self, path: str, pH_min: float = 7.0, pH_max: float = 7.0) -> Chem.Mol:
        """
        Function for loading and protonating molecule given by path to PDB file.
        Achieved with dimorphite_dl.
        """

        print("---- [protonation] Loading molecule and preparing SMILES... ----")

        # loading molecule, removing H and generating SMILES
        CCD_mol = [x for x in Chem.SDMolSupplier(path)][0]
        CCD_mol = Chem.RemoveAllHs(CCD_mol)
        CCD_mol_smiles = Chem.MolToSmiles(CCD_mol)
        print(f"---- [protonation] Base SMILES: {CCD_mol_smiles} ----")

        # protonating SMILES with dimorphite
        prot_list = protonate_smiles(CCD_mol_smiles,
                                    ph_min=pH_min,
                                    ph_max=pH_max,
                                    precision=0.5)
        if not prot_list:
            raise RuntimeError(
                "Dimorphite-DL failed to protonate SMILES")

        # protonated SMILES (choosing from variants, taking first variant)
        protonated_smiles = prot_list[0]
        print(
            f"---- [protonation] Protonated SMILES (pH {pH_min, pH_max}): {protonated_smiles} ----")

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

        return prot_mol_no_h
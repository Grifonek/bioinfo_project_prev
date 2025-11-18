from pathlib import Path

from moleculekit.molecule import Molecule
from moleculekit.tools.preparation import systemPrepare

class ProteinService:
    """
    Service responsible for loading and adding hydrogens to proteins using moleculekit.
    """

    def __init__(self, output_dir: str = "./tmp") -> None:
        self.output_dir =Path(output_dir)
    
    def load_and_hydrate(self, path: str, pH: float = 7.0) -> Path:
        """
        Function for loading and adding hydrogens to protein given by path to PDB file.
        Achieved with moleculekit.
        """

        print("---- [protein] Loading protein... ----")

        # loading protein
        protein = Molecule(path)

        print(f"---- [protein] Adding hydrogens to protein (pH {pH})... ----")

        # adding hydrogens
        protein_prepped, report = systemPrepare(
            protein, return_details=True, pH=pH)

        # creating output PDB file
        out_pdb = self.output_dir / "protein_prepared.pdb"
        protein_prepped.write(str(out_pdb))

        print(f"---- [protein] Protein saved to {out_pdb} ----")
        return out_pdb
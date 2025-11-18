import os
from pathlib import Path

from services.protein_services import ProteinService
from services.protonation_services import ProtonationService
from services.hydride_services import HydrideService

class MainRunner:
    def __init__(self, output_dir: str = "./tmp") -> None:
        self.output_dir = Path(output_dir)
        self.protein_services = ProteinService(output_dir=str(self.output_dir))
        self.protonation_services = ProtonationService()
        self.hydride_services = HydrideService(output_path=str(self.output_dir / "ligand_ph7_h_hydride.pdb"))
    
    def run(self, protein_path: str, ligand_path: str, pH: float = 7.0, pH_min: float = 6.8, pH_max: float = 7.9) -> None:
        # making sure that ./tmp dir exists
        os.makedirs("./tmp", exist_ok=True)

        # 1) loading and hydrating protein
        self.protein_services.load_and_hydrate(protein_path, pH=pH)

        # 2) protonating ligand molecule
        protonated_mol = self.protonation_services.protonate_with_dimorphite(ligand_path, pH_min=pH_min, pH_max=pH_max)
        print(protonated_mol)

        # 3) adding hydrogens to ligand molecule
        self.hydride_services.add_hs_with_hydride(protonated_mol)

if __name__ == "__main__":
    runner = MainRunner(output_dir="./tmp")
    runner.run(
        protein_path="./PDBs/AF-L8BU87-F1-model_v6.pdb",
        ligand_path="./PDBs/PPI_ideal.sdf",
        pH=7.0,
        pH_min=6.8,
        pH_max=7.9
    )
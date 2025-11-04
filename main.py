import os

from services.protein_services import load_and_hydrate_protein
from services.protonation_services import protonate_with_dimorphite
from services.hydride_services import prepare_ligand_from_protonated_rdkit_mol


def main():
    os.makedirs("./tmp", exist_ok=True)

    # 1) loading and hydrating protein
    load_and_hydrate_protein(path="./PDBs/AF-L8BU87-F1-model_v6.pdb", pH=7.0)

    # 2) protonating ligand molecule
    protonated_mol = protonate_with_dimorphite(
        path="./PDBs/BUA.pdb", pH_min=7.0, pH_max=7.0)

    # 3) adding hydrogens to ligand molecule
    prepare_ligand_from_protonated_rdkit_mol(protonated_mol)

    print("Completed successfully. Output saved in ./tmp")
    return


if __name__ == "__main__":
    main()

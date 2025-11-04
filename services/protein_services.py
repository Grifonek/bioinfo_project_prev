from moleculekit.molecule import Molecule
from moleculekit.tools.preparation import systemPrepare


def load_and_hydrate_protein(path: str, pH: float):
    """
    Function for loading and adding hydrogens to protein given by path to PDB file.
    Achieved with moleculekit.
    """

    print("---- Loading protein... ----")

    # loading protein
    protein = Molecule(path)

    print("---- Adding hydrogens to protein (pH 7)... ----")

    # adding hydrogens
    protein_prepped, report = systemPrepare(
        protein, return_details=True, pH=pH)

    # creating output PDB file
    protein_prepped.write("./tmp/protein_prepared.pdb")

    print("---- Protein saved to ./tmp/protein_prepared.pdb ----")
    return

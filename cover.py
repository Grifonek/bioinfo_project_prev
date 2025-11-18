import os
import random

from typing import Tuple

from Bio import PDB
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB import PDBIO

class AtomSelector(PDB.Select):
    """
    Support class for Biopython.
    After initialization, a set with all full ids of the atoms to be written into the substructure must be stored in self.full_ids.
    """
    def accept_atom(self, atom):
        return int(atom.full_id in self.full_ids)

def cover(orig_protein: str, radius: int | float, shift_coords: Tuple[float, float, float]):
    """
    Function for selecting some atom from protein PDB file, selecting substructure within radius around that atom,
    shifting substructure by some Angströms and saving modified protein as PDB file.

    args:
    orig_protein = path to protein PDB file
    radius = radius to search around selected atom
    shift_coords = x,y,z coordinates for moving substructure
    """

    # init
    parser = PDBParser(QUIET=True)
    io = PDBIO()

    # loading protein PDB with biopython
    orig_protein_structure = parser.get_structure("structure", orig_protein)

    # selecting all atoms from original protein structure
    atoms = list(orig_protein_structure.get_atoms())

    # creating kdtree of all atoms for fast searching
    kdtree = NeighborSearch(atoms)

    # creating coords for selected point = atom
    # picking random atom, should be selectable by user later
    center_atom = random.choice(atoms)
    coords = center_atom.coord

    # getting atoms of substructure around atom chosen by coords
    orig_protein_subatoms = kdtree.search(center=coords, radius=radius, level="A")

    # selector for substructure
    selector = AtomSelector()
    selector.full_ids = set([atom.full_id for atom in orig_protein_subatoms])

    # creating folder if does not exist
    os.makedirs("./tmp/cover", exist_ok=True)

    # saving substructure as PDB file
    io.set_structure(orig_protein_structure)
    io.save(file=f"./tmp/cover/substructure.pdb", select=selector, preserve_atom_numbering=True)

    # loading PDB substructure and shift each atom by some Angströms
    substructure = parser.get_structure("substructure", "./tmp/cover/substructure.pdb")
    num_of_shifted_atoms = 0
    for atom in substructure.get_atoms():
        num_of_shifted_atoms += 1
        atom.coord = atom.coord + shift_coords

    # saving shifted substructure in PDB file
    io.set_structure(substructure)
    io.save(file=f"./tmp/cover/shifted_substructure.pdb")

    # loading shifted substructure PDB
    shifted_substructure = parser.get_structure("structure", "./tmp/cover/shifted_substructure.pdb")

    # getting shifted atoms
    shifted_atoms = {atom.full_id: atom.coord for atom in shifted_substructure.get_atoms()}
    matched_atoms = 0
    for atom in orig_protein_structure.get_atoms():
        if atom.full_id in shifted_atoms:
            matched_atoms += 1
            atom.coord = shifted_atoms[atom.full_id]

    # saving modified protein
    io.set_structure(orig_protein_structure)
    io.save(file=f"./tmp/cover/modified_protein.pdb")

    # validation of result
    if matched_atoms == num_of_shifted_atoms:
        print("---- [cover] Substructure shifted and merged successfully. ----")
    else:
        raise ValueError("Something went wrong, number of shifted atoms do not match with number of atoms in substructure!")

    return

cover("./PDBs/AF-L8BU87-F1-model_v6.pdb", 6.0, [1.0, 1.0, 1.0])
#/usr/bin/env pypy
from math import sqrt
from atoms import Atom
from bonds import Bond

class Molecule:
    def __init__(self, molecule):
        self.name = molecule[0]
        self.number_of_atoms = len(molecule[1]["atoms"])
        atoms = []
        for atom in molecule[1]["atoms"]:
            atoms.append(Atom(atom))
        self.atoms = atoms
        bonds = []
        for bond in molecule[1]["bonds"]:
            bonds.append(Bond(bond))
        self.bonds = bonds
        self.total_charge = molecule[1]["total_charge"]
        list_with_bonded_atoms = []
        for bond in self.bonds:
            list_with_bonded_atoms.append(tuple(bond.bonded_atoms))
        self.set_with_bonded_atoms = set(list_with_bonded_atoms)


    @property
    def name(self):
        return self.name

    @property
    def number_of_atoms(self):
        return self.number_of_atoms


    def get_atom_type_with_idx(self, index):
        if index > self.number_of_atoms:
            print("Index is higher than number of atoms.")
            return False
        return self.atoms[index - 1].symbol

    @property
    def formal_charge(self):
        return self.total_charge


    def get_distance_between_atoms(self, index1, index2):
        if index1 > self.number_of_atoms:
            print("Index of atom 1 is higher than number of atoms.")
            return False
        if index2 > self. number_of_atoms:
            print("Index of atom 2 is higher than number of atoms")
            return False
        x1, y1, z1 = self.atoms[index1 - 1].position
        x2, y2, z2 = self.atoms[index2 - 1].position
        distance = sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
        return distance

    def get_bond_type_between_atoms(self, index1, index2):
        if index2 < index1:
            index1, index2 = index2, index1
        if (index1, index2) in self.set_with_bonded_atoms:
            for bond in self.bonds:
                if bond.bonded_atoms == [index1, index2]:
                    return bond.bond_type

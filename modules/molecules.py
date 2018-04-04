from numpy import sqrt, array
from .atoms import Atom
from numba import jit
from scipy import spatial
import numpy as np


class Molecule:
    def __init__(self, molecule, parameters_keys=None):
        self._name = molecule[0]
        self._number_of_atoms = len(molecule[1]["atoms"])
        self._atoms = []
        self._atoms_types = []
        for atom in molecule[1]["atoms"]:
            self._atoms.append(Atom(atom))
            self._atoms_types.append(atom[0])
        self._bonds = []
        self._bonds_types = []
        for bond in molecule[1]["bonds"]:
            self._bonds.append((bond[0], bond[1]))
            self._bonds_types.append(bond[2])
        self._total_charge = molecule[1]["total_charge"]
        self._bonded_atoms = [[] for _ in range(self._number_of_atoms+1)]
        list_of_all_bonds = [[] for _ in range(self._number_of_atoms+1)]
        for x, bond in enumerate(self._bonds):
            atom1 = bond[0]
            atom2 = bond[1]
            self._bonded_atoms[atom1].append(atom2)
            self._bonded_atoms[atom2].append(atom1)
            list_of_all_bonds[atom1].append(self._bonds_types[x])
            list_of_all_bonds[atom2].append(self._bonds_types[x])
        self._highest_bond_of_atoms = []
        for x in list_of_all_bonds[1:]:
            try:
                self._highest_bond_of_atoms.append(max(x))
            except:
                self._highest_bond_of_atoms.append(1)
        self._highest_bond_of_atoms.insert(0, 0)
        self._symbol_gravity = []
        for x, atom in enumerate(self._atoms):
            self._symbol_gravity.append(str(atom.symbol) + "~" + str(self._highest_bond_of_atoms[x+1]))
        self._bonded_bonded_atoms = [[] for _ in range(self._number_of_atoms+1)]
        for index, x in enumerate(self._bonded_atoms[1:]):
            for y in x:
                self._bonded_bonded_atoms[index+1].extend(self._bonded_atoms[y])
            self._bonded_bonded_atoms[index+1] = list(set(self._bonded_bonded_atoms[index+1]))
            try:
                self._bonded_bonded_atoms[index+1].remove(index+1)
            except:
                pass
        atom_cords = array([atom.position for atom in self._atoms])
        self._distance_matrix = spatial.distance.cdist(atom_cords, atom_cords)
        if parameters_keys:
            self._atoms_gravity_bonded_atoms = []
            for atom in self._symbol_gravity:
                self._atoms_gravity_bonded_atoms.append("{}~x".format(atom))
            for index, atom in enumerate(self._atoms):
                for parameter in parameters_keys:
                    if atom.symbol == parameter[0]:
                        searched_atom = parameter[2]
                        for bonded_atom in self.bonded_atoms[index + 1]:
                            if self._atoms_types[bonded_atom-1] == searched_atom:
                                self._atoms_gravity_bonded_atoms[index] = "{}~{}".format(self._symbol_gravity[index], searched_atom)

    def symbol_gravity(self, index):
        return self._symbol_gravity[index-1]

    def symbol_gravity_bonded_atoms(self, index):
        return self._atom_gravity_bonded_atoms[index-1]

    def symbol_to_number(self, atomic_types, type):
        s_numbers = []
        if type == "atom":
            for atom in self._atoms_types:
                s_numbers.append(atomic_types.index(atom))
            self.s_numbers = s_numbers
        elif type == "atom~high_bond":
            for atom_gravity in self._symbol_gravity:
                s_numbers.append(atomic_types.index(atom_gravity))
            self.s_numbers = s_numbers
        elif type == "atom~high_bond~bonded_atoms":
            for atom_gravity_bonded_atoms in self._atom_gravity_bonded_atoms:
                s_numbers.append(atomic_types.index(atom_gravity_bonded_atoms))
            self.s_numbers = s_numbers

    def symbols(self, type):
        if type == "atom":
            return self._atoms_types
        elif type == "atom~high_bond":
            return self._symbol_gravity
        elif type == "atom~high_bond~bonded_atoms":
            return self._atom_gravity_bonded_atoms

    @property
    def symbol_grav(self):
        return self._symbol_gravity

    @property
    def highest_bond_of_atoms(self):
        return self._highest_bond_of_atoms

    @property
    def name(self):
        return self._name

    def __len__(self):
        return self._number_of_atoms

    @property
    def matrix_of_distance(self):
        return self._distance_matrix

    @property
    def atoms_types(self):
        return self._atoms_types

    @property
    def formal_charge(self):
        return self._total_charge

    @property
    def bonded_atoms(self):
        return self._bonded_atoms

    @property
    def bonded_bonded_atoms(self):
        return self._bonded_bonded_atoms

    def set_length_correction(self, correction):
        if correction == 1:
            pass
        else:
            self._distance_matrix *= correction

    def get_atom_type_with_idx(self, index):
        return self._atoms_types[index - 1]


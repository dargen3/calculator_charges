from numpy import sqrt, array
from .atoms import Atom
from numba import jit
from scipy import spatial
import numpy as np

@jit(nopython=True)
def bonded_atoms(index, set_of_bonds):
    bonded_atoms_list = []
    for bond in set_of_bonds:
        if bond[0] == index:
            bonded_atoms_list.append(bond[1])
        elif bond[1] == index:
            bonded_atoms_list.append(bond[0])
    return bonded_atoms_list


@jit(nopython=True)
def distance(x1, y1, z1, x2, y2, z2):
    return sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)


@jit(nopython=True)
def bond_type(index1, index2, bonds):
    if index2 < index1:
        index1, index2 = index2, index1
    for x, bond in enumerate(bonds):
        if bond[0] == index1 and bond[1] == index2:
            return x


class Molecule:
    def __init__(self, molecule):
        self._name = molecule[0]
        self._number_of_atoms = len(molecule[1]["atoms"])
        atoms = []
        atoms_types = []
        for atom in molecule[1]["atoms"]:
            atoms.append(Atom(atom))
            atoms_types.append(atom[0])
        self._atoms_types = atoms_types
        self._atoms = atoms
        bonds = []
        bonds_types = []
        for bond in molecule[1]["bonds"]:
            bonds.append((bond[0], bond[1]))
            bonds_types.append(bond[2])
        self._bonds = bonds
        self._bonds_types = bonds_types
        self._total_charge = molecule[1]["total_charge"]
        list_of_atoms = [[] for _ in range(self._number_of_atoms+1)]
        list_of_all_bonds = [[] for _ in range(self._number_of_atoms+1)]
        for x, bond in enumerate(self._bonds):
            atom1 = bond[0]
            atom2 = bond[1]
            list_of_atoms[atom1].append(atom2)
            list_of_atoms[atom2].append(atom1)
            list_of_all_bonds[atom1].append(self._bonds_types[x])
            list_of_all_bonds[atom2].append(self._bonds_types[x])
        self._bonded_atoms = list_of_atoms
        length = len(sorted(self._bonded_atoms[1:], key=len, reverse=True)[0])
        self._corected_bonded_atoms = np.array([xi+[-1]*(length-len(xi)) for xi in self._bonded_atoms[1:]])
        list_of_highest_bond = []
        for x in list_of_all_bonds[1:]:
            list_of_highest_bond.append(max(x))
        list_of_highest_bond.insert(0, 0)
        self._highest_bond_of_atoms = list_of_highest_bond
        symbol_gravity = []
        for x, atom in enumerate(self._atoms):
            symbol_gravity.append(str(atom.symbol) + "~" + str(self._highest_bond_of_atoms[x+1]))
        self._symbol_gravity = symbol_gravity
        list_of_bonded_atoms = [[] for _ in range(self._number_of_atoms+1)]
        for index, x in enumerate(self._bonded_atoms[1:]):
            for y in x:
                list_of_bonded_atoms[index+1].extend(self._bonded_atoms[y])
            list_of_bonded_atoms[index+1] = list(set(list_of_bonded_atoms[index+1]))
            list_of_bonded_atoms[index+1].remove(index+1)
        self._bonded_bonded_atoms = list_of_bonded_atoms
        length = len(sorted(self._bonded_bonded_atoms[1:], key=len, reverse=True)[0])
        self._corected_bonded_bonded_atoms = np.array([xi+[-1]*(length-len(xi)) for xi in self._bonded_bonded_atoms[1:]])
        atom_cords = array([atom.position for atom in self._atoms])
        self._distance_matrix = spatial.distance.cdist(atom_cords, atom_cords)

    @property
    def c_bonded_atoms(self):
        return self._corected_bonded_atoms


    @property
    def c_bonded_bonded_atoms(self):
        return self._corected_bonded_bonded_atoms


    def symbol_gravity(self, index):
        return self._symbol_gravity[index-1]


    def contain_atom(self, atom_type):
        if atom_type in self._atoms_types:
            return True
        return False


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

    def symbols(self, type):
        if type == "atom":
            return self._atoms_types
        elif type == "atom~high_bond":
            return self._symbol_gravity

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
        return self._atoms[index - 1].symbol

    def get_distance_between_atoms(self, index1, index2):
        x1, y1, z1 = self._atoms[index1 - 1].position
        x2, y2, z2 = self._atoms[index2 - 1].position
        return distance(x1, y1, z1, x2, y2, z2)
    """
    def get_bond_type_between_atoms(self, index1, index2):
        try:
            return self._bonds_types[bond_type(index1, index2, self._bonds)]
        except TypeError:
            pass
    """
    def get_bonded_atoms(self, index):
        return bonded_atoms(index, self._bonds)


from numpy import sqrt, array
from .atoms import Atom
from numba import jit
from scipy import spatial


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
        self.atoms = atoms
        bonds = []
        bonds_types = []
        for bond in molecule[1]["bonds"]:
            bonds.append((bond[0], bond[1]))
            bonds_types.append(bond[2])
        self.bonds = bonds
        self.bonds_types = bonds_types
        self._total_charge = molecule[1]["total_charge"]
        list_of_atoms = [[] for _ in range(self._number_of_atoms+1)]
        list_of_all_bonds = [[] for _ in range(self._number_of_atoms+1)]
        for x, bond in enumerate(self.bonds):
            atom1 = bond[0]
            atom2 = bond[1]
            list_of_atoms[atom1].append(atom2)
            list_of_atoms[atom2].append(atom1)
            list_of_all_bonds[atom1].append(self.bonds_types[x])
            list_of_all_bonds[atom2].append(self.bonds_types[x])
        list_of_highest_bond = []
        for x in list_of_all_bonds[1:]:
            list_of_highest_bond.append(max(x))
        list_of_highest_bond.insert(0, 0)
        self._highest_bond_of_atoms = list_of_highest_bond
        self._bonded_atoms = list_of_atoms
        atom_cords = array([atom.position for atom in self.atoms])
        self._distance_matrix = spatial.distance.cdist(atom_cords, atom_cords)

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

    def set_length_correction(self, correction):
        if correction == 1:
            pass
        else:
            self._distance_matrix *= correction

    def get_atom_type_with_idx(self, index):
        return self.atoms[index - 1].symbol

    def get_distance_between_atoms(self, index1, index2):
        x1, y1, z1 = self.atoms[index1 - 1].position
        x2, y2, z2 = self.atoms[index2 - 1].position
        return distance(x1, y1, z1, x2, y2, z2)

    def get_bond_type_between_atoms(self, index1, index2):
        try:
            return self.bonds_types[bond_type(index1, index2, self.bonds)]
        except TypeError:
            pass

    def get_bonded_atoms(self, index):
        return bonded_atoms(index, self.bonds)


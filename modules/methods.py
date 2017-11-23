import numpy as np
from sys import exit
from termcolor import colored
from math import erf
import warnings
from numba import jit, int64, float64

warnings.filterwarnings("ignore")


class Arciclass:
    @property
    def method_in_parameters(self):
        return self._parameters_method

    @property
    def right_charges_for_parametrization(self):
        return self._list_with_right_charges

    @property
    def length_correction(self):
        return self._length_correction

    def load_parameters(self, file):
        key_words = 0
        with open(file, "r") as file_with_parameters:
            for line in file_with_parameters:
                try:
                    if line.split()[0] in ["<<global>>", "<<key>>", "<<value_symbol>>", "<<value>>", "<<end>>"]:
                        key_words += 1
                except IndexError:
                    pass
        if key_words != 5:
            exit(colored("ERROR! File with parameters is wrong! File not contain keywords <<global>>, <<key>>, "
                         "<<value_symbol>>, <<value>> and <<end>>!", "red"))
        parameters = {}
        with open(file, "r") as file_with_parameters:
            first_line = file_with_parameters.readline().split()
            method = first_line[1]
            second_line = file_with_parameters.readline().split()
            if second_line[0] != "length_type:" or not second_line[1] in ["Angstrom", "Bohr_radius"]:
                exit(colored("ERROR! File with parameters is wrong!\n length_type: Angstrom \n length_type: Bohr_radius"
                             " \n expected!", "red"))
            if second_line[1] == "Angstrom":
                length_correction = 1.0
            if second_line[1] == "Bohr_radius":
                length_correction = 1.8897261245
            third_line = file_with_parameters.readline().split()
            if third_line[0] != "<<global>>":
                exit(colored("ERROR! File with parameters is on third line wrong! <<global>> expected!", "red"))
            global_line = file_with_parameters.readline().split()
            while global_line[0] != "<<key>>":
                parameters[str(global_line[0])] = float(global_line[1])
                global_line = file_with_parameters.readline().split()
            list_with_keys = []
            key = file_with_parameters.readline().split()
            while key[0] != "<<value_symbol>>":
                list_with_keys.append(key[0])
                key = file_with_parameters.readline().split()
            list_with_value_symbol = []
            value_symbol = file_with_parameters.readline().split()
            while value_symbol[0] != "<<value>>":
                list_with_value_symbol.append(value_symbol[0])
                value_symbol = file_with_parameters.readline().split()
            parameters_atoms = []
            value = file_with_parameters.readline().split()
            parameters_keys = []
            while value[0] != "<<end>>":
                parameters_atoms.append(value[0])
                for x in range(len(list_with_keys), len(list_with_value_symbol) + len(list_with_keys)):
                    key_in_parameters = ""
                    for y in range(len(list_with_keys)):
                        key_in_parameters = key_in_parameters + "~" + value[y]
                    parameters_keys.append(key_in_parameters[1:])
                    key_in_parameters = key_in_parameters + "~" + list_with_value_symbol[x - len(list_with_keys)]
                    parameters[key_in_parameters[1:]] = float(value[x])
                value = file_with_parameters.readline().split()
            self.parameters_keys = sorted(list(set(parameters_keys)))
            pattern = ""
            for x in list_with_keys:
                pattern = pattern + "~" + x
            symbol = ""
            for x in list_with_value_symbol:
                symbol = symbol + "/" + x
            pattern = pattern[1:] + "~(" + symbol[1:] + ")"
            parameters_type = ""
            for x in list_with_keys:
                parameters_type = parameters_type + "~" + x
            self._parameters_method, self.parameters_atoms, self.parameters_type, self.parameters_pattern, \
            self._length_correction, self.parameters = method, parameters_atoms, parameters_type[1:], pattern, \
                                                       length_correction, parameters

    @property
    def set_sorted_parameters(self):
        list_of_parameters = []
        list_of_keys = []
        for key in sorted(self.parameters):
            list_of_parameters.append(self.parameters[key])
            list_of_keys.append(key)
        self._sorted_parameters_values = list_of_parameters
        self.sorted_parameters_keys = list_of_keys
        return self.sorted_parameters_keys

    @property
    def sorted_parameters_values(self):
        return self._sorted_parameters_values

    def load_parameters_from_list(self, list_p):
        self.parameters = dict(zip(self.sorted_parameters_keys, list_p))

    def load_charges_for_par(self, file, list_with_atomic_types, set_of_molecules):
        with open(file, "r") as right_charges_file:
            list_with_right_charges = []
            for line in right_charges_file:
                ls = line.split()
                if len(ls) not in {0, 1, 3}:
                    exit(colored("File " + file + " with right charges is wrong! \n", "red"))
                if len(ls) == 3:
                    list_with_right_charges.append(float(ls[2]))
        self._list_with_right_charges = list_with_right_charges
        all_atomic_types = []
        for molecule in set_of_molecules:
            all_atomic_types.extend(molecule.symbols(self.parameters_type))
        self._all_atomic_types = all_atomic_types
        dict_with_right_charges_by_atom_type = {}
        atom_types_in_set = set(self._all_atomic_types)
        for atom in atom_types_in_set:
            dict_with_right_charges_by_atom_type[str(atom)] = []
        for index, charge in enumerate(self._list_with_right_charges):
            dict_with_right_charges_by_atom_type[self._all_atomic_types[index]].append(charge)

        self.dict_with_right_charges_by_atom_type = dict_with_right_charges_by_atom_type
        self.atom_types_in_set = atom_types_in_set

    def set_atomic_types(self, atomic_types):
        self._atomic_types = atomic_types

    @property
    def atomic_types(self):
        return self._atomic_types

    @property
    def all_atomic_types(self):
        return self._all_atomic_types

    def right_charges_for_parameterization_by_atom_types(self, atom_type):
        # try:
        return self.dict_with_right_charges_by_atom_type[str(atom_type)]
        # except:
        #    return False

    def get_key_in_parameters_gravity(self, i, molecule):
        return molecule.symbol_gravity(i)

    def get_key_in_parameters_atom(self, i, molecule):
        return molecule.get_atom_type_with_idx(i)

    def set_parameters_type(self):
        if self.parameters_type == "atom~high_bond":
            self.symbol = self.get_key_in_parameters_gravity
        elif self.parameters_type == "atom":
            self.symbol = self.get_key_in_parameters_atom

    def get_parameter(self, key):
        return self.parameters[key]

    def make_list_of_lists_of_parameters(self):
        l_of_l = [[] for i in range(len(self.parameters_keys))]
        for x in sorted(self.parameters):
            if x[0].isupper():
                for y in self.parameters_keys:
                    if y == x[:len(y)]:
                        l_of_l[self.parameters_keys.index(y)].append(self.parameters[x])
        self.list_of_lists_of_parameters = np.array(l_of_l, dtype=np.float64)


@jit(nopython=True, nogil=True, cache=True)
def eem_calculate(num_of_atoms, kappa, matrix_of_distance, parameters_values, parameters_keys, formal_charge):
    matrix = np.empty((num_of_atoms + 1, num_of_atoms + 1), dtype=np.float64)
    vector = np.empty(num_of_atoms + 1, dtype=np.float64)
    matrix[:num_of_atoms, :num_of_atoms] = kappa / matrix_of_distance
    matrix[num_of_atoms, :] = 1.0
    matrix[:, num_of_atoms] = 1.0
    matrix[num_of_atoms, num_of_atoms] = 0.0
    for i in range(num_of_atoms):
        symbol = parameters_keys[i]
        matrix[i][i] = parameters_values[symbol][1]
        vector[i] = -parameters_values[symbol][0]
    vector[-1] = formal_charge
    results = np.linalg.solve(matrix, vector)
    return results[:-1]


class EEM(Arciclass):
    def calculate(self, molecule):
        return eem_calculate(len(molecule), self.get_parameter("kappa"), molecule.matrix_of_distance,
                             self.list_of_lists_of_parameters, molecule.s_numbers, molecule.formal_charge)

    def calculate_slow(self, molecule):
        num_of_atoms = len(molecule)
        matrix = np.empty((num_of_atoms + 1, num_of_atoms + 1), dtype=float)
        vector = np.empty(shape=[num_of_atoms + 1], dtype=float)
        matrix[:num_of_atoms, :num_of_atoms] = self.get_parameter("kappa") / molecule.matrix_of_distance
        matrix[num_of_atoms, :] = 1.0
        matrix[:, num_of_atoms] = 1.0
        matrix[num_of_atoms, num_of_atoms] = 0.0
        for i in range(1, num_of_atoms + 1):
            symbol = self.symbol(i, molecule)
            i_min1 = i - 1
            matrix[i_min1][i_min1] = self.get_parameter(symbol + "~beta")
            vector[i_min1] = -self.get_parameter(symbol + "~alfa")
        vector[-1] = molecule.formal_charge
        results = np.linalg.solve(matrix, vector)
        return results[:-1]


@jit(nopython=True, nogil=True, cache=True)
def sfkeem_calculate(num_of_atoms, sigma, parameters_values, parameters_keys, matrix_of_distance, formal_charge):
    matrix = np.ones((num_of_atoms + 1, num_of_atoms + 1), dtype=np.float64)
    vector = np.empty(num_of_atoms + 1, dtype=np.float64)
    for x in range(num_of_atoms):
        symbol = parameters_keys[x]
        value = parameters_values[symbol][1]
        matrix[x, :-1] = matrix[x, :-1] * value
        matrix[:-1, x] = matrix[:-1, x] * value
        vector[x] = - parameters_values[symbol][0]
    matrix[:num_of_atoms, :num_of_atoms] = 2.0 * np.sqrt(matrix[:num_of_atoms, :num_of_atoms]) * \
                                           (1 / np.cosh(matrix_of_distance * sigma))
    vector[-1] = formal_charge
    matrix[num_of_atoms, num_of_atoms] = 0.0
    results = np.linalg.solve(matrix, vector)
    return results[:-1]


class SFKEEM(Arciclass):
    def calculate(self, molecule):
        return sfkeem_calculate(len(molecule), self.get_parameter("sigma"), self.list_of_lists_of_parameters,
                                molecule.s_numbers, molecule.matrix_of_distance, molecule.formal_charge)

    def calculate_slow(self, molecule):
        num_of_atoms = len(molecule)
        matrix = np.ones((num_of_atoms + 1, num_of_atoms + 1), dtype=np.float64)
        vector = np.empty(num_of_atoms + 1, dtype=np.float64)
        for x in range(num_of_atoms):
            symbol = self.symbol(x + 1, molecule)
            value = self.get_parameter(symbol + "~beta")
            matrix[x, :-1] = matrix[x, :-1] * value
            matrix[:-1, x] = matrix[:-1, x] * value
            vector[x] = - self.get_parameter(symbol + "~alfa")
        sigma = self.get_parameter("sigma")
        matrix[:num_of_atoms, :num_of_atoms] = 2.0 * np.sqrt(matrix[:num_of_atoms, :num_of_atoms]) * (
            1 / np.cosh(molecule.matrix_of_distance * sigma))
        vector[-1] = molecule.formal_charge
        matrix[num_of_atoms, num_of_atoms] = 0.0
        results = np.linalg.solve(matrix, vector)
        return results[:-1]


@jit(nopython=True)
def coulomb_integral(cor, rad_1, rad_2, distance):
    return cor * erf(np.sqrt(rad_1 * rad_2 / (rad_1 + rad_2)) * distance) / distance


@jit(nopython=True, nogil=True, cache=True)
def qeq_calculate(num_of_atoms, matrix_of_distance, parameters_keys, parameters_values, correlation, formal_charge):
    matrix = np.empty((num_of_atoms + 1, num_of_atoms + 1), dtype=np.float64)
    vector = np.empty(num_of_atoms + 1, dtype=np.float64)
    matrix[:num_of_atoms, :num_of_atoms] = matrix_of_distance
    matrix[num_of_atoms, :] = 1.0
    matrix[:, num_of_atoms] = 1.0
    matrix[num_of_atoms, num_of_atoms] = 0.0
    vector_rad = np.empty(num_of_atoms, dtype=np.float64)
    for x in range(num_of_atoms):
        vector_rad[x] = parameters_values[parameters_keys[x]][2]
    for i in range(num_of_atoms):
        symbol = parameters_keys[i]
        matrix[i][i] = parameters_values[symbol][1]
        vector[i] = -parameters_values[symbol][0]
        for j in range(i + 1, num_of_atoms):
            rad1 = vector_rad[i]
            rad2 = vector_rad[j]
            distance = matrix_of_distance[i][j]
            matrix[i][j] = matrix[j][i] = correlation * erf(np.sqrt(rad1 * rad2 / (rad1 + rad2)) * distance) / distance
    vector[-1] = formal_charge
    results = np.linalg.solve(matrix, vector)
    return results[:-1]


class QEq(Arciclass):
    def calculate(self, molecule):
        return qeq_calculate(len(molecule), molecule.matrix_of_distance, molecule.s_numbers,
                             self.list_of_lists_of_parameters, self.get_parameter("correlation"),
                             molecule.formal_charge)

    def calculate_slow(self, molecule):
        num_of_atoms = len(molecule)
        matrix = np.empty((num_of_atoms + 1, num_of_atoms + 1), dtype=float)
        vector = np.empty(shape=[num_of_atoms + 1], dtype=float)
        matrix[:num_of_atoms, :num_of_atoms] = molecule.matrix_of_distance
        matrix[num_of_atoms, :] = 1.0
        matrix[:, num_of_atoms] = 1.0
        matrix[num_of_atoms, num_of_atoms] = 0.0
        correlation = self.get_parameter("correlation")
        rad = []
        for x in range(1, num_of_atoms + 1):
            rad.append(self.get_parameter(self.symbol(x, molecule) + "~radius"))
        for i in range(num_of_atoms):
            symbol = self.symbol(i + 1, molecule)
            matrix[i][i] = self.get_parameter(symbol + "~hardness")
            vector[i] = -self.get_parameter(symbol + "~electronegativity")
            for j in range(i + 1, num_of_atoms):
                matrix[i][j] = matrix[j][i] = coulomb_integral(correlation, rad[i], rad[j], matrix[i][j])
        vector[-1] = molecule.formal_charge
        results = np.linalg.solve(matrix, vector)
        return results[:-1]


def final(results):
    try:
        end_core = sum(results) / len(results)
        for x in range(len(results)):
            results[x] -= end_core
        return results
    except ZeroDivisionError:
        return results


@jit(nopython=True, nogil=True, cache=True)
def saefm_calculate(num_of_atoms, parameters_values, parameters_keys, bonded_atoms, bonded_bonded_atoms):
    results = np.array([0.0 for x in range(num_of_atoms)])
    for i in range(1, num_of_atoms + 1):
        bonded_atoms_elektronegativity = 0
        symbol_i = parameters_keys[i - 1]
        for atom in bonded_atoms[i - 1]:
            if atom == -1:
                break
            symbol_atom = parameters_keys[atom - 1]
            bonded_atoms_elektronegativity += parameters_values[symbol_atom][0]
        for aatom in bonded_bonded_atoms[i - 1]:
            if aatom == -1:
                break
            symbol_aatom = parameters_keys[aatom - 1]
            bonded_atoms_elektronegativity += parameters_values[symbol_aatom][1]
        charge = bonded_atoms_elektronegativity * parameters_values[symbol_i][3] - parameters_values[symbol_i][2]
        results[i - 1] = charge
    re = np.array([0.0 for x in range(num_of_atoms)])
    for i in range(1, num_of_atoms + 1):
        rec = 0
        for atom in bonded_atoms[i - 1]:
            if atom == -1:
                break
            symbol_atom = parameters_keys[atom - 1]
            rec += parameters_values[symbol_atom][4] * results[atom - 1]
        re[i - 1] = rec + results[i - 1]
    end_core = np.sum(re) / len(re)
    for x in range(len(re)):
        re[x] -= end_core
    return re


class SAEFM(Arciclass):
    def calculate(self, molecule):
        num_of_atom = len(molecule)
        return saefm_calculate(num_of_atom, self.list_of_lists_of_parameters, molecule.s_numbers,
                               molecule.c_bonded_atoms, molecule.c_bonded_bonded_atoms)

    def calculate_(self, molecule):
        num_of_atoms = len(molecule)
        results = []
        for i in range(1, num_of_atoms + 1):
            symbol_i = self.symbol(i, molecule)
            bonded_atoms_elektronegativity = 0
            bonded_atoms = molecule.bonded_atoms[i]
            for atom in bonded_atoms:
                symbol_atom = self.symbol(atom, molecule)
                bonded_atoms_elektronegativity += self.get_parameter(symbol_atom + "~1")
            bonded_bonded_atoms = molecule.bonded_bonded_atoms[i]
            for aatom in bonded_bonded_atoms:
                symbol_aatom = self.symbol(aatom, molecule)
                bonded_atoms_elektronegativity += self.get_parameter(symbol_aatom + "~2")

            charge = bonded_atoms_elektronegativity * self.get_parameter(symbol_i + "~4") - self.get_parameter(
                symbol_i + "~3")
            results.append(charge)
        re = []
        for i in range(1, num_of_atoms + 1):
            rec = 0
            bonded_atoms = molecule.bonded_atoms[i]
            for atom in bonded_atoms:
                symbol_atom = self.symbol(atom, molecule)
                rec += self.get_parameter(symbol_atom + "~5") * results[atom - 1]
            re.append(rec + results[i - 1])
        return final(re)

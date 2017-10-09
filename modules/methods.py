import numpy as np
from sys import exit
from termcolor import colored
from scipy.special import erf
import warnings
from numba import jit
warnings.filterwarnings("ignore")


class arciclass:
    @property
    def _method_in_parameters(self):
        return self._parameters_method

    @property
    def _get_parameterized_atom_type(self):
        try:
            return self._parametrized_atom_type
        except:
            return False

    @property
    def _right_charges_for_parametrization(self):
        return self._list_with_right_charges

    @property
    def _length_correction(self):
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
            exit(colored("ERROR! File with parameters is wrong! File not contain keywords <<global>>, <<key>>, <<value_symbol>>, <<value>> and <<end>>!",               "red"))
        parameters = {}
        with open(file, "r") as file_with_parameters:
            first_line = file_with_parameters.readline().split()
            method = first_line[1]
            second_line = file_with_parameters.readline().split()
            if second_line[0] != "length_type:" or not second_line[1] in ["Angstrom", "Bohr_radius"]:
                exit(colored("ERROR! File with parameters is wrong!\n length_type: Angstrom \n length_type: Bohr_radius \n expected!",                    "red"))
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
            while value[0] != "<<end>>":
                parameters_atoms.append(value[0])
                for x in range(len(list_with_keys), len(list_with_value_symbol) + len(list_with_keys)):
                    key_in_parameters = ""
                    for y in range(len(list_with_keys)):
                        key_in_parameters = key_in_parameters + "~" + value[y]
                    key_in_parameters = key_in_parameters + "~" + list_with_value_symbol[x - len(list_with_keys)]
                    parameters[key_in_parameters[1:]] = float(value[x])
                value = file_with_parameters.readline().split()
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
            parameters_type[1:]
            self._parameters_method, self.parameters_atoms, self.parameters_type, self.parameters_pattern, self._length_correction, self.parameters = method, parameters_atoms, parameters_type[1:], pattern, length_correction, parameters

    def set_sorted_parameters(self):
        list_of_parameters = []
        list_of_keys = []
        for key in sorted(self.parameters):
            list_of_parameters.append(self.parameters[key])
            list_of_keys.append(key)
        self.parameters = list_of_parameters
        self.sorted_parameters_keys = list_of_keys
        self.global_sorted_parameters_keys = list_of_keys
        return self.sorted_parameters_keys

    def set_sorted_parameters_keys_for_atom_type(self, keys):
        self.sorted_parameters_keys = keys

    def set_global_sorted_parameters_keys(self):
        self.sorted_parameters_keys = self.global_sorted_parameters_keys

    def load_parameters_from_list(self, list):
        self.parameters = list

    def load_charges_for_par(self, file):
        with open(file, "r") as right_charges_file:
            list_with_right_charges = []
            for line in right_charges_file:
                l = line.split()
                if len(l) not in {0, 1, 3}:
                    exit(colored("File " + file + " with right charges is wrong! \n", "red"))
                if len(l) == 3:
                    list_with_right_charges.append(float(l[2]))
        self._list_with_right_charges = list_with_right_charges

    def load_charges_for_par_with_molecules(self, file):
        with open(file, "r") as right_charges_file:
            list_with_right_charges = []
            molecules_charges = []
            for line in right_charges_file:
                l = line.split()
                if len(l) not in {0, 1, 3}:
                    exit(colored("File " + file + " with right charges is wrong! \n", "red"))
                if len(l) == 3:
                    molecules_charges.append(float(l[2]))
                if len(l) == 0:
                    list_with_right_charges.append(molecules_charges)
                    molecules_charges = []
            if list_with_right_charges == [] or set(molecules_charges) != list_with_right_charges[-1]:
                list_with_right_charges.append(molecules_charges)
        self._list_with_right_charges_by_molecule = list_with_right_charges

    def load_charges_for_par_by_atom_types(self, file, list_with_atomic_types):
        with open(file, "r") as right_charges_file:
            dict_with_right_charges_by_atom_type = {}
            for atom in list_with_atomic_types:
                dict_with_right_charges_by_atom_type[str(atom)] = []
            num_of_atoms = 0
            for line in right_charges_file:
                l = line.split()
                if len(l) not in {0, 1, 3}:
                    exit(colored("File " + file + " with right charges is wrong! \n", "red"))
                if len(l) == 3:
                    try:
                        dict_with_right_charges_by_atom_type[l[1]].append(float(l[2]))
                    except KeyError:
                        pass
                    num_of_atoms = num_of_atoms + 1
        self.dict_with_right_charges_by_atom_type = dict_with_right_charges_by_atom_type
        return num_of_atoms

    def right_charges_for_parameterization_by_atom_types(self, atom_type):
        try:
            return self.dict_with_right_charges_by_atom_type[str(atom_type)]
        except:
            return False

    def set_parameterized_atom_type(self, atom_type):
        self._parametrized_atom_type = atom_type

    def get_key_in_parameters_gravity(self, i, molecule):
        return molecule.get_atom_type_with_idx(i) + "~" + str(molecule._highest_bond_of_atoms[i])

    def get_key_in_parameters_atom(self, i, molecule):
        return molecule.get_atom_type_with_idx(i)

    def set_parameters_type(self):
        if self.parameters_type == "atom~high_bond":
            self.symbol = self.get_key_in_parameters_gravity
        elif self.parameters_type == "atom":
            self.symbol = self.get_key_in_parameters_atom

    def get_parameter_calc(self, key):
        return self.parameters[key]

    def get_parameter_para(self, key):
        return self.parameters[self.sorted_parameters_keys.index(key)]

    def set_getting_parameters(self, decision):
        if decision == "parameterization":
            self.get_parameter = self.get_parameter_para
        if decision == "calculation":
            self.get_parameter = self.get_parameter_calc


class EEM(arciclass):
    def calculate(self, molecule):
        num_of_atoms = molecule._number_of_atoms
        matrix = np.empty((num_of_atoms + 1, num_of_atoms + 1), dtype=float)
        vector = np.empty(shape=[num_of_atoms + 1], dtype=float)
        matrix[:num_of_atoms, :num_of_atoms] = np.true_divide(np.array([self.get_parameter("kappa")]), molecule._matrix_of_distance)
        matrix[num_of_atoms, :] = 1.0
        matrix[:, num_of_atoms] = 1.0
        matrix[num_of_atoms, num_of_atoms] = 0.0
        for i in range(1, num_of_atoms + 1):
            symbol = self.symbol(i, molecule)
            i_min1 = i-1
            matrix[i_min1][i_min1] = self.get_parameter(symbol + "~beta")
            vector[i_min1] = -self.get_parameter(symbol + "~alfa")
        vector[-1] = molecule._formal_charge
        results = np.linalg.solve(matrix, vector)
        return results[:-1]


@jit(nopython=True)
def sfkeem_formula(dist, p1, p2, sigma):
    return 2.0 * np.sqrt(p1 * p2) * (1.0 / np.cosh(dist * sigma))


class SFKEEM(arciclass):
    def calculate(self, molecule):
        num_of_atoms = molecule._number_of_atoms
        matrix = np.empty((num_of_atoms + 1, num_of_atoms + 1), dtype=float)
        vector = np.empty(shape=[num_of_atoms + 1], dtype=float)
        matrix[:num_of_atoms, :num_of_atoms] = molecule._matrix_of_distance
        matrix[num_of_atoms, :] = 1.0
        matrix[:, num_of_atoms] = 1.0
        matrix[num_of_atoms, num_of_atoms] = 0.0
        list_with_parameters = []
        for i in range(1, num_of_atoms+1):
            symbol = self.symbol(i, molecule)
            list_with_parameters.append(self.get_parameter(symbol + "~beta"))
            vector[i-1] = -(self.get_parameter(symbol + "~alfa"))
        sigma = self.get_parameter("sigma")
        for i in range(num_of_atoms):
            matrix[i][i] = list_with_parameters[i]*2
            for j in range(i+1, num_of_atoms):
                matrix[i][j] = matrix[j][i] = sfkeem_formula(matrix[i][j], list_with_parameters[i], list_with_parameters[j], sigma)
        vector[-1] = molecule._formal_charge
        results = np.linalg.solve(matrix, vector)
        return results[:-1]


class QEq(arciclass):
    def calculate(self, molecule):
        num_of_atoms = molecule._number_of_atoms
        matrix = np.empty((num_of_atoms + 1, num_of_atoms + 1), dtype=float)
        vector = np.empty(shape=[num_of_atoms + 1], dtype=float)
        matrix[:num_of_atoms, :num_of_atoms] = molecule._matrix_of_distance
        matrix[num_of_atoms, :] = 1.0
        matrix[:, num_of_atoms] = 1.0
        matrix[num_of_atoms, num_of_atoms] = 0.0
        correction = self.get_parameter("cor")
        for i in range(num_of_atoms):
            matrix[i][i] = self.get_parameter(self.symbol(i + 1, molecule) + "~hardness")
            vector[i] = -self.get_parameter(self.symbol(i + 1, molecule) + "~electronegativity")
            for j in range(i + 1, num_of_atoms):
                distance = matrix[i][j]
                rad_1 = self.get_parameter(self.symbol(i, molecule) + "~radius")
                rad_2 = self.get_parameter(self.symbol(i, molecule) + "~radius")
                matrix[i][j] = matrix[j][i] = correction * erf(np.sqrt(rad_1 * rad_2 / (rad_1 + rad_2)) * distance) / distance
        vector[-1] = molecule._formal_charge
        results = np.linalg.solve(matrix, vector)
        return results[:-1]


def final(results):
    try:
        end_core = sum(results) / len(results)
        for x in range(len(results)):
            results[x] = results[x] - end_core
        return results
    except ZeroDivisionError:
        return results


class SAEFM(arciclass):
    def calculate(self, molecule, atomic_type):
        num_of_atoms = molecule._number_of_atoms
        results = []
        for i in range(1, num_of_atoms + 1):
            symbol = self.symbol(i, molecule)
            if type(atomic_type) == str and symbol == atomic_type or type(atomic_type) == list and symbol in atomic_type:
                bonded_atoms_elektronegativity = []
                bonded_atoms = molecule._bonded_atoms[i]
                for atom in bonded_atoms:
                    bonded_atoms_elektronegativity.append(self.get_parameter(symbol + "~" + self.symbol(atom, molecule)))
                charge = (sum(bonded_atoms_elektronegativity) - self.get_parameter(symbol + "~" + symbol)) / len(bonded_atoms_elektronegativity)
                results.append(charge)
        if type(atomic_type) == list:
            return final(results)
        return results

class CSAEFM(arciclass):
    def calculate(self, molecule):
        num_of_atoms = molecule._number_of_atoms
        results = []
        for i in range(1, num_of_atoms + 1):
            symbol = self.symbol(i, molecule)
            bonded_atoms = molecule._bonded_atoms[i]
            for atom in bonded_atoms:
                diff = (charge_i - charge_atom)*(parameter1 - parameter2)
            """

            napsat molecule.right_charges
            results = molecule._right_charges


            propojit load_chargeds_by molecule s modulem molecule
            """
        return results




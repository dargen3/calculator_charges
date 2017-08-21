import numpy as np
from sys import exit
from termcolor import colored
from math import erf
import warnings
warnings.filterwarnings("ignore")


class arciclass:
    def symbol_gravity(self, i, num_of_atoms, molecule):
        symbol = molecule.get_atom_type_with_idx(i + 1)
        gravity = 0
        for j in range(num_of_atoms):
            try:
                bond = molecule.get_bond_type_between_atoms(i + 1, j + 1)
                if str(bond) == "SINGLE" and gravity < 1:
                    gravity = 1
                if str(bond) == "DOUBLE" and gravity < 2:
                    gravity = 2
                if str(bond) == "TRIPLE" and gravity < 3:
                    gravity = 3
            except AttributeError:
                continue
        key_in_parameters = symbol + "~" + str(gravity)
        return key_in_parameters

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
            exit(colored(
                "ERROR! File with parameters is wrong! File not contain keywords <<global>>, <<key>>, <<value_symbol>>, <<value>> and <<end>>!",
                "red"))
        parameters = {}
        with open(file, "r") as file_with_parameters:
            first_line = file_with_parameters.readline().split()
            method = first_line[1]
            second_line = file_with_parameters.readline().split()
            if second_line[0] != "length_type:" or not second_line[1] in ["Angstrom", "Bohr_radius"]:
                exit(colored(
                    "ERROR! File with parameters is wrong!\n length_type: Angstrom \n length_type: Bohr_radius \n expected!",
                    "red"))
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
            self.parameters_method, self.parameters_atoms, self.parameters_type, self.parameters_pattern, self.length_correction, self.parameters = method, parameters_atoms, parameters_type[1:], pattern, length_correction, parameters

    def get_parameter(self, key):
        if isinstance(self.parameters, dict):
            return self.parameters[key]
        if isinstance(self.parameters, np.ndarray):
            return self.parameters[self.sorted_parameters_keys.index(key)]

    @property
    def method_in_parameters(self):
        return self.parameters_method

    def set_sorted_parameters(self):
        list_of_parameters = []
        list_of_keys = []
        for key in sorted(self.parameters):
            list_of_parameters.append(self.parameters[key])
            list_of_keys.append(key)
        self.parameters = list_of_parameters
        self.sorted_parameters_keys = list_of_keys

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
        self.list_with_right_charges = list_with_right_charges



    @property
    def right_charges_for_parametrization(self):
        return self.list_with_right_charges


class EEM(arciclass):
    def calculate(self, molecule):
        num_of_atoms = molecule.number_of_atoms
        matrix = np.empty(shape=[num_of_atoms + 1, num_of_atoms + 1], dtype=float)
        vector = np.empty(shape=[num_of_atoms + 1], dtype=float)
        matrix[num_of_atoms] = 1
        matrix[:num_of_atoms] = 1
        matrix[num_of_atoms][num_of_atoms] = 0
        kappa = self.get_parameter("kappa")
        for i in range(num_of_atoms):
            for j in range(i + 1, num_of_atoms):
                distance = molecule.get_distance_between_atoms(i + 1, j + 1) * self.length_correction
                matrix[i][j] = kappa / distance
                matrix[j][i] = kappa / distance
        if self.parameters_type == "atom~high_bond":
            for i in range(num_of_atoms):
                key_in_parameters = self.symbol_gravity(i, num_of_atoms, molecule)
                matrix[i][i] = self.get_parameter(str(key_in_parameters) + "~" + "beta")
                vector[i] = -self.get_parameter(str(key_in_parameters) + "~" + "alfa")
        if self.parameters_type == "atom":
            for i in range(num_of_atoms):
                symbol = molecule.get_atom_type_with_idx(i + 1)
                matrix[i][i] = self.get_parameter(str(symbol) + "~" + "beta")
                vector[i] = -self.get_parameter(str(symbol) + "~" + "alfa")
        vector[-1] = molecule.formal_charge
        results = np.linalg.solve(matrix, vector)
        return results


class SFKEEM(arciclass):
    def calculate(self, molecule):
        num_of_atoms = molecule.number_of_atoms
        matrix = np.empty(shape=[num_of_atoms + 1, num_of_atoms + 1], dtype=float)
        vector = np.empty(shape=[num_of_atoms + 1], dtype=float)
        matrix[num_of_atoms] = 1
        matrix[:num_of_atoms] = 1
        matrix[num_of_atoms][num_of_atoms] = 0
        list_with_parameters = []
        if self.parameters_type == "atom":
            for i in range(num_of_atoms):
                key_in_parameters = molecule.get_atom_type_with_idx(i + 1)
                list_with_parameters.append(self.get_parameter(str(key_in_parameters) + "~" + "beta"))
                vector[i] = -(self.get_parameter(str(key_in_parameters) + "~" + "alfa"))
        if self.parameters_type == "atom~high_bond":
            for i in range(num_of_atoms):
                key_in_parameters = self.symbol_gravity(i, num_of_atoms, molecule)
                matrix[i][i] = self.get_parameter(str(key_in_parameters) + "~" + "beta")
                list_with_parameters.append(self.get_parameter(str(key_in_parameters) + "~" + "beta"))
                vector[i] = -self.get_parameter(str(key_in_parameters) + "~" + "alfa")
        sigma = self.get_parameter("sigma")
        for i in range(num_of_atoms):
            for j in range(num_of_atoms):
                distance = molecule.get_distance_between_atoms(i + 1, j + 1) * self.length_correction
                matrix[i][j] = 2.0 * np.sqrt(list_with_parameters[i] * list_with_parameters[j]) * (1.0 / np.cosh(distance * sigma))
        vector[-1] = molecule.formal_charge
        results = np.linalg.solve(matrix, vector)
        return results

class QEq(arciclass):
    def calculate(self, molecule):
        num_of_atoms = molecule.number_of_atoms
        matrix = np.empty(shape=[num_of_atoms + 1, num_of_atoms + 1], dtype=float)
        vector = np.empty(shape=[num_of_atoms + 1], dtype=float)
        matrix[num_of_atoms] = 1
        matrix[:num_of_atoms] = 1
        matrix[num_of_atoms][num_of_atoms] = 0
        radius_of_atoms = []
        for i in range(num_of_atoms):
            atom = molecule.get_atom_type_with_idx(i + 1)
            radius_of_atoms.append(self.get_parameter(str(atom) + "~radius"))
        max_distance = 2.0 * np.sqrt(-np.log(1.e-9)/min(radius_of_atoms))
        for i in range(num_of_atoms):
            for j in range(i + 1, num_of_atoms):
                distance = molecule.get_distance_between_atoms(i + 1, j + 1) * self.length_correction
                if distance > max_distance:
                    matrix[i][j] = 1.0 / distance
                    matrix[j][i] = 1.0 / distance
                else:
                    rad_1 = radius_of_atoms[i]
                    rad_2 = radius_of_atoms[j]
                    #print rad_1,rad_2
                    value = erf(np.sqrt(rad_1 * rad_2 / (rad_1 + rad_2)) * distance) / distance
                    print value, 1.0/distance
                    matrix[i][j] = value
                    matrix[j][i] = value
        for i in range(num_of_atoms):
            symbol = molecule.get_atom_type_with_idx(i + 1)
            matrix[i][i] = self.get_parameter(str(symbol) + "~" + "hardness")
            vector[i] = self.get_parameter(str(symbol) + "~" + "electronegativity")
        vector[-1] = molecule.formal_charge
        #pprint(matrix)
        #pprint(vector)
        results = np.linalg.solve(matrix, vector)
        return results

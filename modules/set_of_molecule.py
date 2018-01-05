from .molecules import Molecule
from sys import exit
from termcolor import colored
from tabulate import tabulate


class Set_of_molecule:
    def __init__(self, file, parameters_keys=None):
        if parameters_keys:
            parameters_keys = [parameter.split("~") for parameter in parameters_keys if parameter[-1] != "x"]
        with open(file, 'r') as if_sdf:
            num_of_molecules = 0
            num_of_lines = 0
            for line in if_sdf:
                num_of_lines += 1
                if line.strip() == "$$$$":
                    num_of_molecules += 1
            self._num_of_lines = num_of_lines
            if line.strip() != "$$$$":
                exit(colored("ERROR! End of sdf file must be: $$$$", "red"))
        list_with_all_molecules = []
        actual_mol = 0
        with open(file, 'r') as sdf_file:
            while actual_mol < num_of_molecules:
                name = sdf_file.readline().strip()
                sdf_file.readline()
                sdf_file.readline()
                line = sdf_file.readline()
                if str(line[-5]) == "2":
                    number_of_atoms = int(line[:3])
                    number_of_bonds = int(line[3:6])
                    list_of_atoms = []
                    for x in range(number_of_atoms):
                        line = sdf_file.readline().split()
                        list_of_atoms.append((line[3], x + 1, float(line[0]), float(line[1]), float(line[2])))
                    list_of_bond = []
                    for x in range(number_of_bonds):
                        line = sdf_file.readline()
                        type_of_bond = int(line[8])
                        atom1 = int(line[:3])
                        atom2 = int(line[3:6])
                        if atom2 < atom1:
                            atom1, atom2 = atom2, atom1
                        list_of_bond.append((atom1, atom2, type_of_bond))
                    line = sdf_file.readline().strip()
                    charge = 0
                    while line != "$$$$":
                        line = sdf_file.readline().strip()
                        try:
                            ls = line.split()
                            if ls[0] == "M" and ls[1] == "CHG":
                                for x in range(5, 3 + int(ls[2]) * 2, 2):
                                    charge = charge + float(ls[x - 1])
                        except IndexError:
                            pass
                    dict_with_mol_data = {"atoms": list_of_atoms, "bonds": list_of_bond, "total_charge": charge}
                    list_with_all_molecules.append(Molecule((name, dict_with_mol_data), parameters_keys=parameters_keys))
                    actual_mol += 1
                elif str(line[-5]) == "3":
                    list_of_atoms = []
                    for x in range(3):
                        sdf_file.readline()
                    line = sdf_file.readline().split()
                    charge = 0
                    while line[2] != "END" and line[3] != "ATOM":
                        for x in line:
                            if x[:3] == "CHG":
                                charge = charge + float(x[4:])
                        list_of_atoms.append((line[3], int(line[2]), float(line[4]), float(line[5]), float(line[6])))
                        line = sdf_file.readline().split()
                    sdf_file.readline()
                    line = sdf_file.readline().split()
                    list_of_bond = []
                    while line[2] != "END" and line[3] != "BOND":
                        try:
                            type_of_bond = int(line[3])
                        except IndexError:
                            for x in range(10):
                                print(colored(sdf_file.readline(), "red"))
                            exit(colored("ERROR! Wrong type of bond. Multiplicity of bond is " + str(line[2]) +
                                         ". In molecule number: " + str(actual_mol), "red"))
                        atom1 = int(line[4])
                        atom2 = int(line[5])
                        if atom2 < atom1:
                            atom1, atom2 = atom2, atom1
                        list_of_bond.append((atom1, atom2, type_of_bond))
                        line = sdf_file.readline().split()
                    while line != "$$$$":
                        line = sdf_file.readline().strip()
                    actual_mol += 1
                    dict_with_mol_data = {"atoms": list_of_atoms, "bonds": list_of_bond, "total_charge": charge}
                    list_with_all_molecules.append(Molecule((name, dict_with_mol_data)))
                else:
                    for x in range(10):
                        print(colored(sdf_file.readline(), "red"))
                    exit(colored("ERROR! Sdf file is not correct!", "red"))
        self.list_with_molecules = list_with_all_molecules
        self.file = file

    def __len__(self):
        return len(self.list_with_molecules)

    def __getitem__(self, index):
        return self.list_with_molecules[index]

    def __iter__(self):
        return [x for x in self.list_with_molecules]

    @property
    def num_of_lines(self):
        return self._num_of_lines

    def statistics_data(self, write_to_file=None, quiet=False):
        number_of_mol = len(self.list_with_molecules)
        atom_types = []
        for molecule in self.list_with_molecules:
            atom_types.extend(molecule.symbols("atom~high_bond"))
        number_of_atoms = len(atom_types)
        atom_types_counting = {}
        for atom in atom_types:
            if atom not in atom_types_counting:
                atom_types_counting[atom] = 0
            atom_types_counting[atom] += 1
        tab_of_atoms = []
        for atom in sorted(atom_types_counting):
            tab_of_atoms.append((atom, atom_types_counting[atom], round(
                atom_types_counting[atom] / (number_of_atoms / 100), 2)))
        if quiet:
            return set(atom_types)
        if write_to_file == None:
            print("\n\n\nStatistics data from set of molecules from {}\n".format(self.file))
            print("Number of molecules:   {}".format(number_of_mol))
            print("Number of atoms:       {}".format(number_of_atoms))
            print("Number of atoms type:  {}\n".format(len(atom_types_counting)))
            print(tabulate(tab_of_atoms, headers=["Type", "Number", "%"]))
            print("\n\n\n")
        else:
            with open(write_to_file, "w") as file_info:
                file_info.write("\n\n\nStatistics data from set of molecules from {}\n".format(self.file))
                file_info.write("Number of molecules:   {}\n".format(number_of_mol))
                file_info.write("Number of atoms:       {}\n".format(number_of_atoms))
                file_info.write("Number of atoms type:  {}\n\n".format(len(atom_types_counting)))
                file_info.write(tabulate(tab_of_atoms, headers=["Type", "Number", "%"]))
                file_info.write("\n\n\n")
        return set(atom_types)


"""
format of molecule in self.list_with_molecules

('NSC_102026',
 {'atoms': [('C', 1, -0.6392, -0.5604, 0.0011),
            ('C', 2, -0.0547, 1.6192, -0.001),
            ('C', 3, 0.4564, 0.2999, -0.0),
            ('O', 4, 2.1798, -1.2657, 0.0021),
            ('C', 5, 1.8718, -0.0884, 0.0005),
            ('N', 6, -1.7472, 0.2083, 0.0009),
            ('N', 7, -1.3558, 1.5513, -0.0004),
            ('N', 8, 2.8305, 0.8589, -0.0008),
            ('N', 9, -0.5984, -1.9399, -0.0028),
            ('C', 10, -3.1281, -0.2808, 0.0011),
            ('H', 11, 0.5322, 2.5258, -0.0014),
            ('H', 12, 2.5851, 1.7973, -0.002),
            ('H', 13, 3.7659, 0.6022, -0.0004),
            ('H', 14, -1.4238, -2.4495, -0.0058),
            ('H', 15, 0.2556, -2.3999, -0.0026),
            ('H', 16, -3.4706, -0.4017, 1.0288),
            ('H', 17, -3.7677, 0.4367, -0.513),
            ('H', 18, -3.1735, -1.2412, -0.5124)],
  'bonds': [(1, 3, 'DOUBLE'),
            (1, 6, 'SINGLE'),
            (1, 9, 'SINGLE'),
            (2, 3, 'SINGLE'),
            (2, 7, 'DOUBLE'),
            (2, 11, 'SINGLE'),
            (3, 5, 'SINGLE'),
            (4, 5, 'DOUBLE'),
            (5, 8, 'SINGLE'),
            (6, 7, 'SINGLE'),
            (6, 10, 'SINGLE'),
            (8, 12, 'SINGLE'),
            (8, 13, 'SINGLE'),
            (9, 14, 'SINGLE'),
            (9, 15, 'SINGLE'),
            (10, 16, 'SINGLE'),
            (10, 17, 'SINGLE'),
            (10, 18, 'SINGLE')],
  'total_charge': 0})

"""

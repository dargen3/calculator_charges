import os.path
from termcolor import colored
from sys import stdin
from subprocess import call
from os import getcwd, rename, remove
from tempfile import NamedTemporaryFile
from tabulate import tabulate




def control_if_arguments_files_exist_for_par(right_charges, sdf_input, parameters, new_parameters, force, chg_output):
    if not os.path.isfile(right_charges):
        exit(colored("There is no charges file with name " + right_charges + "\n", "red"))
    if not os.path.isfile(sdf_input):
        exit(colored("There is no sdf file with name " + sdf_input + "\n", "red"))
    if not os.path.isfile(parameters):
        exit(colored("There is no parameters file with name " + parameters + "\n", "red"))
    if os.path.isfile(new_parameters):
        if not force:
            print(colored("Warning. There is some file with have the same name like your new parameters file!", "red"))
            print("If you want to replace exist file, please write yes and press enter. Else press enter.")
            decision = stdin.readline().rstrip('\n')
            if decision == "yes":
                print(colored("Exist file will be replaced.\n\n\n", "green"))
            else:
                exit("\n")
    if os.path.isfile(chg_output):
        if not force:
            print(colored("Warning. There is some file with have the same name like your chg output!", "red"))
            print("If you want to replace exist file, please write yes and press enter. Else press enter.")
            decision = stdin.readline().rstrip('\n')
            if decision == "yes":
                os.remove(chg_output)
                print(colored("Exist file was removed.\n\n\n", "green"))
            else:
                print("\n\n")
                exit("\n")
        else:
            os.remove(chg_output)


def control_of_missing_atoms(set_of_molecule, method, parameters):
    list_of_atoms = []
    for molecule in set_of_molecule:
        num_of_atoms = len(molecule)
        for i in range(num_of_atoms):
            symbol = molecule.get_atom_type_with_idx(i)
            if symbol not in list_of_atoms:
                list_of_atoms.append(symbol)
    missing_atoms = set(list_of_atoms) - set(method.parameters_atoms)
    if len(missing_atoms) > 0:
        print(colored("Warning. File with parameters not contain parameters for atom:", "red"))
        for atom in missing_atoms:
            print(atom)
        print("If you want to add missing parameters, please write yes to open parameters file in nano and press enter."
              " Else press enter.")
        decision = stdin.readline().rstrip('\n')
        if decision == "yes":
            call(["nano", parameters])
            method.load_parameters(parameters)
        else:
            exit("\n")
    return list(set(list_of_atoms))


def writing_new_parameters(parameters, new_parameters_file, res, method):
    random_filename = NamedTemporaryFile(prefix="auxiliary_file.", dir=getcwd(), delete=False)
    with open(parameters, "r") as parameters:
        with open(random_filename.name, "w") as new_parameters:
            actual_line = parameters.readline()
            new_parameters.write(actual_line)
            while actual_line.split()[0] != "<<global>>":
                actual_line = parameters.readline()
                new_parameters.write(actual_line)
            actual_line = parameters.readline()
            while actual_line.split()[0] != "<<key>>":
                new_parameters.write(str(actual_line.split()[0]) + " "),
                new_parameters.write(str(float(
                    "{0:.4f}".format(res[method.sorted_parameters_keys.index(actual_line.split()[0])]))) + "\n")
                actual_line = parameters.readline()
            new_parameters.write(actual_line)
            key_list = []
            actual_line = parameters.readline()
            new_parameters.write(actual_line)
            while actual_line.split()[0] != "<<value_symbol>>":
                key_list.append(actual_line.split()[0])
                actual_line = parameters.readline()
                new_parameters.write(actual_line)
            value_symbols_list = []
            actual_line = parameters.readline()
            new_parameters.write(actual_line)
            while actual_line.split()[0] != "<<value>>":
                value_symbols_list.append(actual_line.split()[0])
                actual_line = parameters.readline()
                new_parameters.write(actual_line)
            actual_line = parameters.readline()
            while actual_line.split()[0] != "<<end>>":
                key_of_sorted_parameters = ""
                for x in range(len(key_list)):
                    new_parameters.write(actual_line.split()[x] + "   "),
                    key_of_sorted_parameters = key_of_sorted_parameters + "~" + actual_line.split()[x]
                for x in range(len(key_list), len(key_list) + len(value_symbols_list)):
                    key = key_of_sorted_parameters[1:] + "~" + value_symbols_list[x - len(key_list)]
                    new_parameters.write(str(float("{0:.4f}".format(res[method.sorted_parameters_keys.index(key)]))) +
                                         "   ")
                actual_line = parameters.readline()
                new_parameters.write("\n")
            new_parameters.write("<<end>>\n\n")
    rename(random_filename.name, new_parameters_file)

def write_to_para(parameters, sdf_input, method_par, choised_num_of_mol, table_for_all_atoms, table_for_all_molecules, statistics_data, now, validation):
    with open(parameters, "a") as parameters:
        parameters.write("\n\n\nParameterized set of molecules: " + sdf_input + "\n")
        parameters.write("Date of parameterization: {}\n".format(now.strftime("%Y-%m-%d %H:%M")))
        parameters.write("Method of parameterization: " + method_par)
        parameters.write("\nNumber of parameterized molecules: " + str(choised_num_of_mol) + "\n")
        if validation:
            parameters.write("Mode: Validation 70:30\n\n")
        else:
            parameters.write("Mode: Full set parameterization")
        parameters.write("Statistics for all atoms:\n")
        parameters.write(tabulate(table_for_all_atoms,
                                  headers=["RMSD", "max deviation", "average deviation", "pearson**2",
                                           "num. of atoms"]))
        parameters.write("\n\n")
        parameters.write("Statistics for molecules:\n")
        parameters.write(tabulate(table_for_all_molecules,
                                  headers=["RMSD", "max deviation", "average deviation", "pearson**2",
                                           "num. of molecules"]))
        parameters.write("\n\n")
        parameters.write("Statistics for atomic type:\n")
        parameters.write(tabulate(statistics_data, headers=["atomic type", "RMSD", "max deviation",
                                                            "average deviation", "pearson**2",
                                                            "num. of atoms"]))
        parameters.write("\n\n")



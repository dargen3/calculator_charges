import os.path
from termcolor import colored
from sys import stdin


def writing_to_list(molecule, results):
    num_of_atoms = len(molecule)
    list_with_molecule_data = [molecule.name, str(num_of_atoms)]
    for i in range(num_of_atoms):
        result = results[i]
        list_with_molecule_data.append([i + 1, molecule.get_atom_type_with_idx(i + 1), result])
    list_with_molecule_data.append("")
    return list_with_molecule_data


def writing_to_file(list_with_data, chg_output):
    with open(chg_output, "a") as file_with_results:
        for molecule in list_with_data:
            for item in molecule:
                if isinstance(item, list):
                    file_with_results.write(
                        '{0:>3} {1:>3} {2:>15}'.format(item[0], item[1], str(float("{0:.6f}".format(item[2]))) + "\n"))
                else:
                    file_with_results.write(item + "\n")


def control_if_arguments_files_exist_for_calc(parameters, sdf_input, chg_output, force):
    if not os.path.isfile(parameters):
        exit(colored("There is no parameters file with name " + parameters + "\n", "red"))
    if not os.path.isfile(sdf_input):
        exit(colored("There is no sdf file with name " + sdf_input + "\n", "red"))
    if os.path.isfile(chg_output):
        if force:
            pass
        else:
            print(colored("Warning. There is some file with have the same name like your chg output!", "red"))
            print("If you want to replace exist file, please write yes and press enter. Else press enter.")
            decision = stdin.readline().rstrip('\n')
            if decision == "yes":
                os.remove(chg_output)
                print(colored("Exist file was removed.\n\n\n", "green"))
            else:
                print("\n\n")
                exit(1)

import os.path
from termcolor import colored
from sys import stdin
from modules.set_of_molecule import Set_of_molecule
import importlib


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
        if not force:
            print(colored("Warning. There is some file with have the same name like your chg output!", "red"))
            print("If you want to replace exist file, please write yes and press enter. Else press enter.")
            decision = stdin.readline().rstrip('\n')
            if decision == "yes":
                os.remove(chg_output)
                print(colored("Exist file was removed.\n\n\n", "green"))
            else:
                exit("\n")


def calculate_charges(parameters, sdf_input, chg_output, rewriting_with_force, args_method, logger):
    control_if_arguments_files_exist_for_calc(parameters, sdf_input, chg_output, rewriting_with_force)
    list_with_data = []
    logger.info("Loading molecule data from {} ...".format(sdf_input))
    setm = Set_of_molecule(sdf_input)
    number_of_molecules = len(setm)
    logger.info(colored("Loading molecule data from was successful.".format(sdf_input), "green"))
    logger.info("{} molecules was loaded. \n\n\n".format(number_of_molecules))
    try:
        method = getattr(importlib.import_module("modules.methods"), str(args_method))
    except AttributeError:
        exit(colored("ERROR! Method do not exist or you do not define method!\n\n\n", "red"))
    logger.info("Loading parameters from {} ...".format(parameters))
    method = method()
    method.load_parameters(parameters)
    logger.info("File with parameters: \n---------------------------------------")
    logger.info(colored(open(parameters, 'r').read(), "yellow"))
    logger.info("---------------------------------------")
    if method.method_in_parameters != args_method:
        exit(colored("ERROR! These parameters are for method {} but you want to calculate charges by method {}!" /
                     "\n\n\n".format(method.method_in_parameters, args_method), "red"))
    logger.info(colored("Loading of parameters was sucessfull.\n\n\n", "green"))
    logger.info("Calculating of charges ...")
    method.set_parameters_type()
    num_of_err = 0
    sorted_atomic_types = method.parameters_keys
    for molecule in setm[:number_of_molecules]:
        molecule.symbol_to_number(sorted_atomic_types, method.parameters_type)
    method.make_list_of_lists_of_parameters()
    for molecule in setm[:number_of_molecules]:
        molecule.set_length_correction(method.length_correction)
        try:
            results = method.calculate(molecule)
        except linalg.linalg.LinAlgError:
            print(colored("Molecule {} is not calculated. Solution do not exist.".format(molecule.name), "red"))
            num_of_err += 1
            continue
        except KeyError:
            print(colored("Molecule {} contain atoms, for which are not parameters in file: {}".format(
                molecule.name, parameters), "red"))
            num_of_err += 1
            continue
        list_with_data.append(writing_to_list(molecule, results))
    logger.info(colored("Calculation was successful.", "green"))
    logger.info("{} molecules from {} loaded molecules was calculated. \n\n\n".format(
        number_of_molecules - num_of_err, number_of_molecules))
    logger.info("Writing calculated charges to {} ...".format(chg_output))
    writing_to_file(list_with_data, chg_output)
    logger.info(colored("Writing to {} was successful.\n\n\n".format(chg_output), "green"))

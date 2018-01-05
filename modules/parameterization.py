import os.path
from termcolor import colored
from sys import stdin
from subprocess import call
from os import getcwd, rename, remove
from tempfile import NamedTemporaryFile
from tabulate import tabulate
from .set_of_molecule import Set_of_molecule
import importlib
from scipy.optimize import minimize, differential_evolution, basinhopping
from numba import jit
from scipy import stats
from numpy import linalg, sqrt
from math import isnan
from .calculation import writing_to_list, writing_to_file
from .make_html import make_html
from .alarm import alarm
from .comparison import control_if_arguments_files_exist_for_com, making_dictionary_with_charges_para, \
    statistics_for_all_atoms, plotting, making_final_list, statistics_for_molecules, statistics_for_atom_type
from matplotlib import pyplot as plt
import datetime
from .statistics import statistics as statistics_sdf_and_chg_file


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
    list_of_symbols = []
    for molecule in set_of_molecule:
        symbols = molecule.symbols(method.parameters_type)
        list_of_symbols.extend(symbols)
    missing_atoms = set(list_of_symbols) - set(method.parameters_keys)
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


def write_to_para(parameters, sdf_input, method_par, choised_num_of_mol, table_for_all_atoms, table_for_all_molecules,
                  statistics_data, now, validation):
    with open(parameters, "a") as parameters:
        parameters.write("\n\n\nParameterized set of molecules: " + sdf_input + "\n")
        parameters.write("Date of parameterization: {}\n".format(now.strftime("%Y-%m-%d %H:%M")))
        parameters.write("Method of parameterization: " + method_par)
        parameters.write("\nNumber of parameterized molecules: " + str(choised_num_of_mol) + "\n")
        if validation:
            parameters.write("Mode: Validation 70:30\n\n")
        else:
            parameters.write("Mode: Full set parameterization\n\n")
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


@jit(nopython=True)
def statistics(deviation_list):
    rmsd_list = []
    abs_dev = 0
    for atom_data in deviation_list:
        rmsd_list.append((atom_data[0] - atom_data[1]) ** 2)
        dev = abs(atom_data[0] - atom_data[1])
        if dev > abs_dev:
            abs_dev = dev
    return rmsd_list, abs_dev


def calculating_charges(list_of_parameters, method, set_of_molecule, objective_function):
    method.load_parameters_from_list(list_of_parameters)
    method.make_list_of_lists_of_parameters()
    list_with_results = []
    for molecule in set_of_molecule:
        try:
            list_with_results.extend(method.calculate(molecule))
        except linalg.linalg.LinAlgError:
            pass
    # greater_rmsd = False
    pearson2 = stats.pearsonr(list_with_results, method.right_charges_for_parametrization[:len(list_with_results)])[
                   0] ** 2
    deviation_list = zip(list_with_results, method.right_charges_for_parametrization)
    atomic_types = sorted(method.atom_types_in_set)

    dict_with_right_charges_by_atom_type = {}
    for atom in atomic_types:
        dict_with_right_charges_by_atom_type[str(atom)] = []
    all_atomic_types = method.all_atomic_types

    for x, result in enumerate(list_with_results):
        dict_with_right_charges_by_atom_type[all_atomic_types[x]].append(result)
    statistic_list_rmsd = [0] * len(atomic_types)
    statistic_list_pear = [0] * len(atomic_types)
    statistic_list_absd = [0] * len(atomic_types)
    for atom in atomic_types:
        deviation_list_at = zip(dict_with_right_charges_by_atom_type[atom],
                                method.right_charges_for_parameterization_by_atom_types(atom))
        try:
            rmsd_list_at, abs_dev_at = statistics(list(deviation_list_at))
        except ValueError:
            print(list(deviation_list_at))
            exit(colored("ERROR!!! This is few molecules for parameterization", "red"))
        statistic_list_rmsd[atomic_types.index(atom)] = sqrt((1.0 / len(rmsd_list_at)) * sum(rmsd_list_at))
        statistic_list_pear[atomic_types.index(atom)] = stats.pearsonr(dict_with_right_charges_by_atom_type[atom],
                                                        method.right_charges_for_parameterization_by_atom_types(
                                                        atom)[:len(dict_with_right_charges_by_atom_type[atom])])[0] ** 2
        statistic_list_absd[atomic_types.index(atom)] = abs_dev_at
        greater_rmsd = max(statistic_list_rmsd)

    rmsd_list, abs_dev = statistics(list(deviation_list))
    rmsd = sqrt((1.0 / len(rmsd_list)) * sum(rmsd_list))
    if isnan(rmsd):
        rmsd = 1000.0
    if greater_rmsd and isnan(greater_rmsd):
        greater_rmsd = 1000.0
    if isnan(pearson2):
        pearson2 = 0.0
    print("Actual RMSD: " + str(rmsd)[:8] + "   Actual max deviation: " + str(abs_dev)[:8] + "   Worst RMSD: " +
          str(greater_rmsd)[:8] + "   Worst pearson**2: " + str(min(statistic_list_pear))[:8], end="\r")
    if objective_function == "rmsd":
        return rmsd
    elif objective_function == "all":
        return (1 - min(statistic_list_pear)) + (1 - pearson2) + greater_rmsd * 10 + rmsd * 20 + sum(
            statistic_list_pear) / len(statistic_list_pear) + sum(statistic_list_rmsd) / len(statistic_list_rmsd)


def parameterize(args_method, parameters, sdf_input, num_of_parameterized_mol, validation, right_charges,
                 method_parameterization, new_parameters, chg_output, all_mol_to_log, logger, rewriting_with_force,
                 args_save_fig, args_make_html, alarm_after_para):
    control_if_arguments_files_exist_for_par(right_charges, sdf_input, parameters,
                                             new_parameters, rewriting_with_force, chg_output)

    statistics_sdf_and_chg_file(right_charges, sdf_input, logger)
    try:
        method = getattr(importlib.import_module("modules.methods"), args_method)
    except AttributeError:
        exit(colored("ERROR! Method do not exist or you do not define method!\n", "red"))
    logger.info("Loading parameters from {} ...".format(parameters))
    method = method()
    method.load_parameters(parameters)
    logger.info("File with parameters: \n---------------------------------------")
    logger.info(colored(open(parameters, 'r').read(), "yellow"))
    logger.info("---------------------------------------")
    if method.method_in_parameters != args_method:
        exit(colored("ERROR! These parameters are for method {} but you want to calculate charges by method {}!\n\n\n"
                     .format(method.method_in_parameters, args_method), "red"))
    method.set_parameters_type()
    logger.info(colored("Loading of parameters was sucessfull.\n\n\n", "green"))
    logger.info("Loading molecule data from {} ...".format(sdf_input))
    if method.parameters_type == "atom~high_bond~bonded_atoms":
        setm = Set_of_molecule(sdf_input, parameters_keys=method.parameters_keys)
    else:
        setm = Set_of_molecule(sdf_input)
    number_of_molecules = len(setm)
    set_of_molecule = setm[:number_of_molecules]
    all_atomic_types = []
    for molecule in set_of_molecule:
        all_atomic_types.extend(molecule.symbols(method.parameters_type))
    if num_of_parameterized_mol:
        choised_num_of_mol = num_of_parameterized_mol
    else:
        choised_num_of_mol = number_of_molecules
    if validation:
        high_limit = choised_num_of_mol
        choised_num_of_mol = round(choised_num_of_mol * 0.7)
    set_of_molecule = setm[:choised_num_of_mol]
    control_of_missing_atoms(set_of_molecule, method, parameters)
    sorted_atomic_types = method.parameters_keys
    method.set_atomic_types(sorted_atomic_types)
    logger.info(colored("Loading molecule data from {} was successful.\n\n\n".format(sdf_input),
                        "green"))
    logger.info("Loading charges data from {} ...".format(right_charges))
    method.load_charges_for_par(right_charges, setm[:number_of_molecules])
    logger.info(colored("Loading of charges data was sucessfull. \n\n\n", "green"))
    logger.info(str(number_of_molecules) + " molecules was loaded.\n\n\n")
    input_parameters_list = method.sorted_parameters_values
    if args_method == "EEM":
        bounds = [(0.0001, 2)] * len(input_parameters_list)
    elif args_method == "SAEFM":
        bounds = [(-2, 2)] * len(input_parameters_list)
    elif args_method == "QEq":
        bounds = [(0.0001, 4)] * len(input_parameters_list)
    elif args_method == "SFKEEM":
        bounds = [(0.0001, 4)] * len(input_parameters_list)
    else:
        bounds = [(-4, 4)] * len(input_parameters_list)
    for molecule in set_of_molecule:
        molecule.symbol_to_number(sorted_atomic_types, method.parameters_type)
        molecule.set_length_correction(method.length_correction)
    method.make_list_of_lists_of_parameters()
    print("Parameterization running for {} molecules ...\n".format(choised_num_of_mol))
    if method_parameterization == "diferential_evolution":
        res = differential_evolution(calculating_charges, bounds, maxiter=5, popsize=50, polish=False, args=(method, set_of_molecule, "all"))
    elif method_parameterization == "basinhopping":
        res = basinhopping(calculating_charges, input_parameters_list, minimizer_kwargs={"args":(method, set_of_molecule, "all"), "method":"SLSQP"}, stepsize=0.1)
    elif method_parameterization == "combined":
        res = differential_evolution(calculating_charges, bounds, maxiter=5, popsize=50, polish=False, args=(method, set_of_molecule, "rmsd"))
        fun1 = res.fun
        res = minimize(calculating_charges, res.x, method="SLSQP", bounds=bounds,
                       args=(method, set_of_molecule, "rmsd"))
        fun2 = res.fun
        while abs(fun1 - fun2) > fun2 / 1000.0:
            fun1 = fun2
            res = minimize(calculating_charges, res.x, method="SLSQP", bounds=bounds,
                           args=(method, set_of_molecule, "rmsd"))
            fun2 = res.fun
        res = minimize(calculating_charges, res.x, method="SLSQP", bounds=bounds,
                       args=(method, set_of_molecule, "all"))
        fun1 = fun2
        fun2 = res.fun
        while abs(fun1 - fun2) > fun2 / 1000.0:
            fun1 = fun2
            res = minimize(calculating_charges, res.x, method="SLSQP", bounds=bounds,
                           args=(method, set_of_molecule, "all"))
            fun2 = res.fun
    else:
        res = minimize(calculating_charges, input_parameters_list, method="SLSQP", bounds=bounds,
                       args=(method, set_of_molecule, "all"))
    final_parameters = res.x
    writing_new_parameters(parameters, new_parameters, final_parameters, method)
    logger.info("\n\n\nFile with new parameters: \n---------------------------------------")
    logger.info(colored(open(new_parameters, 'r').read(), "yellow"))
    logger.info("---------------------------------------\n\n")
    if validation:
        set_of_molecule = setm[choised_num_of_mol:high_limit]
        for molecule in set_of_molecule:
            molecule.symbol_to_number(sorted_atomic_types, method.parameters_type)
            molecule.set_length_correction(method.length_correction)
    else:
        set_of_molecule = setm[:number_of_molecules]
        if choised_num_of_mol != number_of_molecules:
            for molecule in set_of_molecule[choised_num_of_mol:]:
                molecule.symbol_to_number(sorted_atomic_types, method.parameters_type)
                molecule.set_length_correction(method.length_correction)
    method.make_list_of_lists_of_parameters()
    list_with_data = []
    all_atomic_types_calc = []
    for molecule in set_of_molecule:
        results = method.calculate(molecule)
        list_with_data.append(writing_to_list(molecule, results))
        all_atomic_types_calc.extend(molecule.symbols(method.parameters_type))
    logger.info("Writing calculated charges to {} ...".format(chg_output))
    writing_to_file(list_with_data, chg_output)
    logger.info(colored("Writing to {} was successful.\n\n\n".format(chg_output), "green"))
    plt.switch_backend('agg')
    fig = plt.figure(figsize=(11, 9))
    fig_all = fig.add_subplot(111)
    save_fig = chg_output[:-4] + ".png"
    charges = chg_output
    control_if_arguments_files_exist_for_com(charges, right_charges, save_fig, args_save_fig,
                                             all_mol_to_log, rewriting_with_force)
    if all_mol_to_log:
        with open(all_mol_to_log, "a") as mol_log:
            mol_log.write("name   rmsd   max.dev.   av.dev   pearson**2\n")
    char2 = making_dictionary_with_charges_para(charges, all_atomic_types_calc)
    char1 = making_dictionary_with_charges_para(right_charges, all_atomic_types)
    list_with_atomic_data, dict_with_atomic_data, dictionary_with_molecular_data, atoms = making_final_list(char1,
                                                                                                            char2)
    print("\nStatistics for all atoms:")
    table_for_all_atoms, axis_range, number_of_atoms = statistics_for_all_atoms(list_with_atomic_data,
                                                                                fitting=True)
    RMSD = table_for_all_atoms[0][0]
    pearson_2 = table_for_all_atoms[0][3]
    print(tabulate(table_for_all_atoms, headers=["RMSD", "max deviation", "average deviation", "pearson**2",
                                                 "num. of atoms"]))
    print("\n\nStatistics for molecules:")
    table_for_all_molecules = statistics_for_molecules(dictionary_with_molecular_data, all_mol_to_log)
    print(tabulate(table_for_all_molecules, headers=["RMSD", "max deviation", "average deviation", "pearson**2",
                                                     "num. of molecules"]))
    statistics_data = []
    atoms = sorted(atoms)
    for atom in atoms:
        statistics_data.append(statistics_for_atom_type(atom, dict_with_atomic_data[atom], atoms,
                                                        fig_all, charges, right_charges,
                                                        args_save_fig, axis_range, save_fig[:-4]))
    print("\n\nStatistics for atomic type:")
    if method_parameterization:
        method_par = method_parameterization
    else:
        method_par = "minimize"
    now = datetime.datetime.now()
    write_to_para(new_parameters, sdf_input, method_par, choised_num_of_mol, table_for_all_atoms,
                  table_for_all_molecules, statistics_data, now, validation)
    print(tabulate(statistics_data, headers=["atomic type", "RMSD", "max deviation", "average deviation",
                                             "pearson**2", "num. of atoms"]))
    print("\n\n\n")
    if args_make_html:
        make_html(charges[:-4], sdf_input, method, table_for_all_atoms[0],
                  table_for_all_molecules[0], save_fig, atoms, statistics_data,
                  now.strftime("%Y-%m-%d %H:%M"), method_parameterization, choised_num_of_mol,
                  validation)
    if alarm_after_para:
        alarm()
    plotting(charges, right_charges, args_save_fig, save_fig,
             fig_all, fig, RMSD, pearson_2, axis_range, number_of_atoms)

#!/usr/bin/python3
from sys import exit, stdin
import argparse
import argcomplete
from numpy import linalg, sqrt, array, concatenate
import importlib
from termcolor import colored
from subprocess import call, check_call, CalledProcessError
from scipy.optimize import minimize, basinhopping, differential_evolution
from scipy import stats
from matplotlib import pyplot as plt
from tabulate import tabulate
import warnings
from modules.set_of_molecule import Set_of_molecule
import modules.func_calc as calc
import modules.func_para as para
import modules.func_comp as comp
import modules.func_stat as stat
import modules.func_html as make_html
import logging
from math import isnan
from numba import jit
from multiprocessing import Pool
import datetime

def settings_argparse():
    global args
    global logger
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", help="It can be chosen calculation, parameterization, comparison or statistics.",
                        required=True, choices=("calculation", "parameterization", "comparison", "statistics"))
    parser.add_argument("--method", help="Method to calculating charges or for parameterization.")
    parser.add_argument("--sdf_input", help="Sdf file with molecules data.")
    parser.add_argument("--parameters", help="File with parameters.")
    parser.add_argument("--all_mol_to_log", help="For comparison. Results of all molecules are saved into log file.")
    parser.add_argument("--chg_output", help="Output chg file with charges.")
    parser.add_argument("--right_charges", help="File with charges from QM.")
    parser.add_argument("--new_parameters", help="File to save parameters from parameterization.")
    parser.add_argument("--charges", help="File with charges.")
    parser.add_argument("--num_of_parameterized_mol", help="Only first N molecule will be parameterized.")
    parser.add_argument("--save_fig", action="store_true", help="Save figures of comparison.")
    parser.add_argument("--validation", action="store_true",
                        help="From set of molecules will be 70% used for parameterization and 30% for validation.")
    parser.add_argument("--choised_molecules", help="PRIVATE ARGUMENT!")
    parser.add_argument("--start_from_ones", action="store_true",
                        help="All parameters are before parameterization set to 1.")
    parser.add_argument("--atom_types_for_para", help="Only for SAEFM!", nargs="+")
    parser.add_argument("--method_parameterization",
                        help="It can be chosen minimize, basinhopping or diferential_evolution.",
                        choices=("minimize", "basinhopping", "diferential_evolution"))
    parser.add_argument("--alarm_after_para", action="store_true",
                        help="Alarm after parameterization. You need instalated bash's play!")
    parser.add_argument("--make_html", action="store_true", help="Make html after comparison. Only for parameterization.")
    parser.add_argument("--statistic_on_end_of_para",
                        help="PRIVATE ARGUMENT! Statistics data from para. is added to end of parameters file.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity.")
    parser.add_argument("-f", "--rewriting_with_force", action="store_true",
                        help="All existed files with the same names like your outputs will be replaced.")
    argcomplete.autocomplete(parser)
    args = parser.parse_args()
    if args.mode == "calculation":
        if args.method is None or args.sdf_input is None or args.parameters is None or args.chg_output is None:
            parser.error("For calculation must be choisen --method, --sdf_input, --parameters and --chg_output!")
    elif args.mode == "parameterization":
        if args.right_charges is None or args.method is None or args.sdf_input is None or args.parameters is None or \
                        args.new_parameters is None or args.chg_output is None:
            parser.error("For parameterization must be choisen --right_charges, --method, --sdf_input, " +
                         "--parameters, --chg_output and --new_parameters!")
    elif args.mode == "comparison":
        if args.charges is None or args.right_charges is None:
            parser.error("For comparison must be choisen --charges and --right_charges!")
    elif args.mode == "statistics":
        if args.charges is None or args.sdf_input is None:
            parser.error("For statistics must be choisen --charges and --sdf_input")
    if args.verbose:
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        logger.addHandler(ch)
    else:
        logger = logging.getLogger()
        logger.setLevel(logging.ERROR)
        ch = logging.StreamHandler()
        ch.setLevel(logging.ERROR)
        logger.addHandler(ch)
    logger.info("\n\n")
    logger.info("---------------------------------------")
    logger.info("\n\n")
    logger.info(colored("Verbosity turned on.\nMode:" + args.mode + "\n\n\n", "blue"))
    if args.verbose is False:
        print("\n\n")


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


def molecule_calculate(molecules):
    """
    if args.method[0:2] == "SA":
        result = method.calculate(molecule, method.get_parameterized_atom_type)
    else:
        result = method.calculate(molecule)
    return result
    """
    list_with_results = []
    for molecule in molecules:
        if args.method[0:2] == "SA":
            result = method.calculate(molecule, method.get_parameterized_atom_type)
        else:
            result = method.calculate(molecule)
        list_with_results.extend(result)
    return list_with_results

def split_list(arr, size):
    arrs = []
    while len(arr) > size:
        pice = arr[:size]
        arrs.append(pice)
        arr = arr[size:]
    arrs.append(arr)
    return arrs


def calculating_charges(list_of_parameters):

    method.load_parameters_from_list(list_of_parameters)
    """
    result = pool.map(molecule_calculate, split_list(set_of_molecule, int(len(set_of_molecule)/2) + 1))
    list_with_results = list(concatenate(result))
    """
    list_with_results = []
    for molecule in set_of_molecule:
        if args.method[0:2] == "SA":
            result = method.calculate(molecule, method.get_parameterized_atom_type)
        else:
            result = method.calculate(molecule)
        list_with_results.extend(result)

    greater_rmsd = False
    if method.get_parameterized_atom_type is False:
        pearson2 = stats.pearsonr(list_with_results, method.right_charges_for_parametrization[:len(list_with_results)])[
                       0] ** 2
        deviation_list = zip(list_with_results, method.right_charges_for_parametrization)
        atomic_types = sorted(method.atomic_types)
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
            rmsd_list_at, abs_dev_at = statistics(list(deviation_list_at))
            statistic_list_rmsd[atomic_types.index(atom)] = sqrt((1.0 / len(rmsd_list_at)) * sum(rmsd_list_at))
            statistic_list_pear[atomic_types.index(atom)] = stats.pearsonr(dict_with_right_charges_by_atom_type[atom],
                                                                           method.right_charges_for_parameterization_by_atom_types(
                                                                               atom)[:len(
                                                                               dict_with_right_charges_by_atom_type[
                                                                                   atom])])[0] ** 2
            statistic_list_absd[atomic_types.index(atom)] = abs_dev_at
            greater_rmsd = max(statistic_list_rmsd)
    else:
        deviation_list = zip(list_with_results, method.right_charges_for_parameterization_by_atom_types(
            method.get_parameterized_atom_type))
        pearson2 = stats.pearsonr(list_with_results,
                                  method.right_charges_for_parameterization_by_atom_types(
                                      method.get_parameterized_atom_type)[:len(list_with_results)])[0] ** 2
    rmsd_list, abs_dev = statistics(list(deviation_list))
    rmsd = sqrt((1.0 / len(rmsd_list)) * sum(rmsd_list))
    if isnan(rmsd):
        rmsd = 1000.0
    if greater_rmsd and isnan(greater_rmsd):
        greater_rmsd = 1000.0
    if isnan(pearson2):
        pearson2 = 1000.0
    if method.get_parameterized_atom_type is False:
        print("Actual RMSD: " + str(rmsd)[:9] + "     Actual max deviation: " + str(abs_dev)[:9] + "     Worst RMSD: " +
              str(greater_rmsd)[:9] + "     Worst pearson**2: " + str(min(statistic_list_pear))[:9], end="\r")
        return rmsd #(1 - min(statistic_list_pear)) + (1 - pearson2) + greater_rmsd * 10 + rmsd * 10 + sum(
        #    statistic_list_pear) / len(statistic_list_pear) + sum(statistic_list_rmsd) / len(statistic_list_rmsd)
    print("Actual RMSD: " + str(rmsd) + "         Actual max deviation: " + str(
        abs_dev) + "         Actual pearson**2: " + str(pearson2), end="\r")
    return rmsd  + abs_dev / 10.0 + (1-pearson2)


if __name__ == "__main__":
    settings_argparse()
    warnings.filterwarnings("ignore")
    if args.mode == "calculation":
        calc.control_if_arguments_files_exist_for_calc(args.parameters, args.sdf_input, args.chg_output,
                                                       args.rewriting_with_force)
        list_with_data = []
        logger.info("Loading molecule data from " + str(args.sdf_input) + " ...")
        setm = Set_of_molecule(args.sdf_input)
        number_of_molecules = len(setm)
        logger.info(colored("Loading molecule data from " + args.sdf_input + " was successful.", "green"))
        logger.info(str(number_of_molecules) + " molecules was loaded. \n\n\n")
        try:
            method = getattr(importlib.import_module("modules.methods"), str(args.method))
        except AttributeError:
            exit(colored("ERROR! Method do not exist or you do not define method!\n\n\n", "red"))
        logger.info("Loading parameters from " + str(args.parameters) + " ...")
        method = method()
        method.load_parameters(args.parameters)
        logger.info("File with parameters: \n---------------------------------------")
        logger.info(colored(open(args.parameters, 'r').read(), "yellow"))
        logger.info("---------------------------------------")
        if method.method_in_parameters != args.method:
            exit(colored("ERROR! These parameters are for method " + str(
                method.method_in_parameters) + " but you want to calculate charges by method " +
                         str(args.method) + "!\n\n\n", "red"))
        logger.info(colored("Loading of parameters was sucessfull.\n\n\n", "green"))
        logger.info("Calculating of charges ...")
        method.set_parameters_type()
        num_of_err = 0
        if args.method == "SAEFM":
            atomic_types = para.control_of_missing_atoms(setm[:number_of_molecules], method, args.parameters)
        if args.validation:
            set_of_molecule = set_of_molecule[int(args.choised_molecules):int(int(args.choised_molecules) * 1 / 0.7)]
        sorted_atomic_types = method.parameters_keys
        for molecule in setm[:number_of_molecules]:
            molecule.symbol_to_number(sorted_atomic_types, method.parameters_type)
        method.make_list_of_lists_of_parameters()
        for molecule in setm[:number_of_molecules]:
            molecule.set_length_correction(method.length_correction)
            try:
                if args.method == "SAEFM":
                    results = method.calculate(molecule, atomic_types)
                else:
                    results = method.calculate(molecule)
            except linalg.linalg.LinAlgError:
                print(colored("Molecule " + molecule.name + " is not calculated. Solution do not exist.", "red"))
                num_of_err += 1
                continue
            except KeyError:
                print(colored("Molecule " + molecule.name + " contain atoms, for which are not parameters in file: " +
                              args.parameters, "red"))
                num_of_err += 1
                continue
            list_with_data.append(calc.writing_to_list(molecule, results))
        logger.info(colored("Calculation was successful.", "green"))
        logger.info(str(number_of_molecules - num_of_err) + " molecules from " + str(
            number_of_molecules) + " loaded molecules was calculated. \n\n\n")
        logger.info("Writing calculated charges to " + args.chg_output + " ...")
        calc.writing_to_file(list_with_data, args.chg_output)
        logger.info(colored("Writing to " + args.chg_output + " was successful.\n\n\n", "green"))

    elif args.mode == "parameterization":
        para.control_if_arguments_files_exist_for_par(args.right_charges, args.sdf_input, args.parameters,
                                                      args.new_parameters, args.rewriting_with_force)
        call_stat = "./calculator_charges.py --mode statistics --sdf_input " + args.sdf_input + \
                    "  --charges " + args.right_charges
        if args.verbose:
            call_stat += " -v "
        exit_value_s = call(call_stat, shell=True)
        if exit_value_s != 0:
            exit()
        logger.info("---------------------------------------\n\n\n")
        try:
            method = getattr(importlib.import_module("modules.methods"), str(args.method))
        except AttributeError:
            exit(colored("ERROR! Method do not exist or you do not define method!\n", "red"))
        logger.info("Loading parameters from " + str(args.parameters) + " ...")
        method = method()
        method.load_parameters(args.parameters)
        logger.info("File with parameters: \n---------------------------------------")
        logger.info(colored(open(args.parameters, 'r').read(), "yellow"))
        logger.info("---------------------------------------")
        if method.method_in_parameters != args.method:
            exit(colored("ERROR! These parameters are for method " + str(
                method.method_in_parameters) + " but you want to calculate charges by method " + str(
                args.method) + "!\n\n\n", "red"))
        method.set_parameters_type()
        logger.info(colored("Loading of parameters was sucessfull.\n\n\n", "green"))
        logger.info("Loading molecule data from " + str(args.sdf_input) + " ...")
        setm = Set_of_molecule(args.sdf_input)
        number_of_molecules = len(setm)
        set_of_molecule = setm[:number_of_molecules]
        for molecule in set_of_molecule:
            molecule.set_length_correction(method.length_correction)
        if args.num_of_parameterized_mol:
            choised_num_of_mol = int(args.num_of_parameterized_mol)
        else:
            choised_num_of_mol = number_of_molecules
        if args.validation:
            choised_num_of_mol = int(choised_num_of_mol * 0.7)
        global pool
        if args.method == "SAEFM":
            set_of_molecule = setm[:choised_num_of_mol]
            a_types = para.control_of_missing_atoms(set_of_molecule, method, args.parameters)
            if args.atom_types_for_para:
                for x in args.atom_types_for_para:
                    if x not in a_types:
                        exit(colored("ERROR! In chooised molecules are no atoms with symbol " + x, "red"))
                atomic_types = args.atom_types_for_para
            else:
                atomic_types = a_types
            logger.info(colored("Loading molecule data from " + str(args.sdf_input) + " was successful.\n\n\n",
                                "green"))
            global_sorted_parameters = method.set_sorted_parameters()
            logger.info("Loading charges data from " + str(args.right_charges) + " ...")
            method.load_charges_for_par_by_atom_types(args.right_charges, atomic_types)
            logger.info(colored("Loading of charges data was sucessfull. \n\n\n", "green"))
            logger.info(str(number_of_molecules) + " molecules was loaded.\n\n\n")
            global_input_parameters_list = method.parameters
            final_parameters = global_input_parameters_list
            logger.info("Parameterization running for " + str(choised_num_of_mol) + " molecules.\n\n\n")
            for atomic_type in atomic_types:
                method.set_parameterized_atom_type(atomic_type)
                input_parameters_list = []
                sorted_parameters = []
                for x in range(len(global_sorted_parameters)):
                    char = ""
                    for character in global_sorted_parameters[x]:
                        if character == "~":
                            break
                        else:
                            char = char + character
                    if char == atomic_type:
                        input_parameters_list.append(global_input_parameters_list[x])
                        sorted_parameters.append(global_sorted_parameters[x])
                method.set_sorted_parameters_keys_for_atom_type(sorted_parameters)
                logger.info("Parameterization running for atom: " + atomic_type)
                bounds = [(-1, 5)] * len(input_parameters_list)
                if args.validation or args.start_from_ones:
                    input_parameters_list = array([1] * len(input_parameters_list))
                if args.method_parameterization == "basinhopping":
                    res = basinhopping(calculating_charges, input_parameters_list, niter=1, T=0.1, stepsize=0.05)
                elif args.method_parameterization == "diferential_evolution":
                    res = differential_evolution(calculating_charges, bounds, maxiter=1)
                else:
                    res = minimize(calculating_charges, input_parameters_list, method="SLSQP", bounds=bounds)
                for x in range(len(res.x)):
                    final_parameters[global_sorted_parameters.index(sorted_parameters[x])] = res.x[x]
                logger.info("\n\n")
            if not args.method_parameterization:
                args.method_parameterization = "minimize"
            method.set_global_sorted_parameters_keys()
        else:
            set_of_molecule = setm[:choised_num_of_mol]
            atomic_types = para.control_of_missing_atoms(set_of_molecule, method, args.parameters)
            method.set_atomic_types(atomic_types)
            logger.info(colored("Loading molecule data from " + str(args.sdf_input) + " was successful.\n\n\n",
                                "green"))
            method.set_sorted_parameters()
            logger.info("Loading charges data from " + str(args.right_charges) + " ...")
            method.load_charges_for_par(args.right_charges)
            method.load_charges_for_par_by_atom_types(args.right_charges, atomic_types)
            logger.info(colored("Loading of charges data was sucessfull. \n\n\n", "green"))
            logger.info(str(number_of_molecules) + " molecules was loaded.\n\n\n")
            input_parameters_list = method.parameters
            if args.validation or args.start_from_ones:
                input_parameters_list = array([1] * len(input_parameters_list))
            bounds = [(-5, 4)] * len(input_parameters_list)
            set_of_molecule = setm[:choised_num_of_mol]
            sorted_atomic_types = method.parameters_keys
            for molecule in set_of_molecule:
                molecule.symbol_to_number(sorted_atomic_types, method.parameters_type)
            method.make_list_of_lists_of_parameters()
            print("Parameterization running for " + str(choised_num_of_mol) + " molecules ...")
            if args.method_parameterization == "basinhopping":
                res = basinhopping(calculating_charges, input_parameters_list, niter=1, T=0.1, stepsize=0.05)
            elif args.method_parameterization == "diferential_evolution":
                res = differential_evolution(calculating_charges, bounds, maxiter=1)
            else:
                res = minimize(calculating_charges, input_parameters_list, method="SLSQP", bounds=bounds)
                args.method_parameterization = "minimize"
            final_parameters = res.x
        para.writing_new_parameters(args.parameters, args.new_parameters, final_parameters, method)
        logger.info("\n\n\nFile with parameters: \n---------------------------------------")
        logger.info(colored(open(args.new_parameters, 'r').read(), "yellow"))
        logger.info("---------------------------------------")
        if args.alarm_after_para:
            try:
                check_call("play modules/alarm/alarm.mp3 2> /dev/null", shell=True)
            except KeyboardInterrupt:
                pass
            except CalledProcessError:
                print(colored("You cannot run mp3 by bash command play!", "red"))
        call1 = "./calculator_charges.py --mode calculation --parameters " + args.new_parameters + \
                " --sdf_input " + args.sdf_input + " --method " + args.method + " --chg_output " + args.chg_output
        call2 = "./calculator_charges.py --mode comparison --charges " + args.chg_output + " --right_charges " + \
                args.right_charges + " --statistic_on_end_of_para " + args.new_parameters + " --sdf_input " + \
                args.sdf_input + " --method_parameterization " + args.method_parameterization + \
                " --choised_molecules " + str(choised_num_of_mol)
        if args.validation:
            call1 += " --validation --choised_molecules " + str(choised_num_of_mol)
        if args.make_html:
            call2 += " --make_html --method " + str(args.method)
        if args.verbose:
            call1 += " -v"
            call2 += " -v"
        if args.rewriting_with_force:
            call1 += " -f"
            call2 += " -f"
        if args.save_fig:
            call2 = call2 + " --save_fig "
        exit_value = call(call1, shell=True)
        if exit_value == 0:
            call(call2, shell=True)

    elif args.mode == "comparison":
        now = datetime.datetime.now()
        fig = plt.figure(figsize=(11, 9))
        fig_all = fig.add_subplot(111)
        save_fig = args.charges[:-4] + ".png"
        comp.control_if_arguments_files_exist_for_com(args.charges, args.right_charges, save_fig,
                                                      args.all_mol_to_log, args.rewriting_with_force)
        if args.all_mol_to_log is not None:
            with open(args.all_mol_to_log, "a") as mol_log:
                mol_log.write("name   rmsd   max.dev.   av.dev   pearson**2\n")
        logger.info("Loading data from " + args.charges + " and " + args.right_charges + " ...")
        char2 = comp.making_dictionary_with_charges(args.charges)
        char1 = comp.making_dictionary_with_charges(args.right_charges)
        logger.info(
            colored("Loading data from " + args.charges + " and " + args.right_charges + " was sucessfull.\n\n\n.",
                    "green"))
        list_with_atomic_data, dictionary_with_molecular_data, atoms = comp.making_final_list(char1, char2)
        print("Statistics for all atoms:")
        table_for_all_atoms, axis_range, number_of_atoms = comp.statistics_for_all_atoms(list_with_atomic_data,
                                                                                         fitting=True)
        RMSD = table_for_all_atoms[0][0]
        pearson_2 = table_for_all_atoms[0][3]
        print(tabulate(table_for_all_atoms, headers=["RMSD", "max deviation", "average deviation", "pearson**2",
                                                     "num. of atoms"]))
        if args.statistic_on_end_of_para:
            with open(args.statistic_on_end_of_para, "a") as parameters:
                parameters.write("\n\n\nParameterized set of molecules: " + args.sdf_input + "\n")
                if args.method_parameterization:
                    method_par = args.method_parameterization
                else:
                    method_par = "minimize"
                parameters.write("Date of parameterization: {}\n".format(now.strftime("%Y-%m-%d %H:%M")))
                parameters.write("Method of parameterization: " + method_par)
                parameters.write("\nNumber of parameterized molecules: " + args.choised_molecules + "\n\n")
                parameters.write("Statistics for all atoms:\n")
                parameters.write(tabulate(table_for_all_atoms,
                                          headers=["RMSD", "max deviation", "average deviation", "pearson**2",
                                                   "num. of atoms"]))
                parameters.write("\n\n")
        print("\n\nStatistics for molecules:")
        table_for_all_molecules = comp.statistics_for_molecules(dictionary_with_molecular_data, args.all_mol_to_log)
        print(tabulate(table_for_all_molecules, headers=["RMSD", "max deviation", "average deviation", "pearson**2",
                                                         "num. of molecules"]))
        if args.statistic_on_end_of_para:
            with open(args.statistic_on_end_of_para, "a") as parameters:
                parameters.write("Statistics for molecules:\n")
                parameters.write(tabulate(table_for_all_molecules,
                                          headers=["RMSD", "max deviation", "average deviation", "pearson**2",
                                                   "num. of molecules"]))
                parameters.write("\n\n")
        statistics_data = []
        atoms = sorted(atoms)
        for atom in atoms:
            statistics_data.append(comp.statistics_for_atom_type(atom, list_with_atomic_data, atoms,
                                                                 fig_all, args.charges, args.right_charges,
                                                                 args.save_fig, axis_range, save_fig[:-4]))
        print("\n\nStatistics for atomic type:")
        print(tabulate(statistics_data, headers=["atomic type", "RMSD", "max deviation", "average deviation",
                                                 "pearson**2", "num. of atoms"]))
        if args.statistic_on_end_of_para:
            with open(args.statistic_on_end_of_para, "a") as parameters:
                parameters.write("Statistics for atomic type:\n")
                parameters.write(tabulate(statistics_data, headers=["atomic type", "RMSD", "max deviation",
                                                                    "average deviation", "pearson**2",
                                                                    "num. of atoms"]))
                parameters.write("\n\n")
        print("\n\n\n")
        comp.plotting(args.charges, args.right_charges, args.save_fig, save_fig,
                      fig_all, fig, RMSD, pearson_2, axis_range, number_of_atoms)
        if args.make_html:
            make_html.make_html(args.charges[:-4], args.sdf_input, args.method, table_for_all_atoms[0], table_for_all_molecules[0], save_fig, atoms, statistics_data, now.strftime("%Y-%m-%d %H:%M"), method_par, args.choised_molecules)
    elif args.mode == "statistics":
        stat.control_if_arguments_files_exist_for_stat(args.charges, args.sdf_input)
        charges_data, sdf_data, charges_names, sdf_names, setm, number_of_lines = stat.is_the_same(args.charges,
                                                                                                   args.sdf_input,
                                                                                                   logger)
        if charges_data == sdf_data and sdf_names == charges_names:
            logger.info(colored("Sdf file is for the same molecules like charge file!", "green"))
        else:
            if charges_data == sdf_data and sdf_names != charges_names:
                print(colored(
                    "Sdf file is for the same molecules like charge file, but molecules is in different order!\n\n\n",
                    "red"))
            elif not (set(sdf_names) & set(charges_names)):
                exit(colored("No together molecules for chg and sdf file.", "red"))
            else:
                print(colored("Sdf file is for diferent set of molecules like charges file!", "red"))
            print(
                "If you want to create new sdf and chg file only with together molecules write yes. Else press enter.")
            decision = stdin.readline().rstrip('\n')
            if decision == "yes":
                stat.create_new_files(sdf_names, charges_names, args.sdf_input, args.charges, setm, number_of_lines)
                exit(colored("Cresting new files was successful!", "green"))
            else:
                exit(1)

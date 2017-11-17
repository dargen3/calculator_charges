#!/usr/bin/python3
from sys import exit, stdin
import seaborn as sns
from numpy import linalg, sqrt, array, concatenate, random
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
import modules.func_clust as clust
import modules.func_html as make_html
import modules.func_argparse as arg
import modules.func_alarm as alarm
from math import isnan
from numba import jit
import datetime
from operator import itemgetter


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


def calculating_charges(list_of_parameters):
    method.load_parameters_from_list(list_of_parameters)
    method.make_list_of_lists_of_parameters()
    list_with_results = []
    for molecule in set_of_molecule:
        try:
            list_with_results.extend(method.calculate(molecule))
        except linalg.linalg.LinAlgError:
            pass
    greater_rmsd = False
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
        try:
            rmsd_list_at, abs_dev_at = statistics(list(deviation_list_at))
        except ValueError:
            exit(colored("ERROR!!! This is few molecules for parameterization", "red"))
        statistic_list_rmsd[atomic_types.index(atom)] = sqrt((1.0 / len(rmsd_list_at)) * sum(rmsd_list_at))
        statistic_list_pear[atomic_types.index(atom)] = stats.pearsonr(dict_with_right_charges_by_atom_type[atom],
                                                                       method.right_charges_for_parameterization_by_atom_types(
                                                                           atom)[:len(
                                                                           dict_with_right_charges_by_atom_type[
                                                                               atom])])[0] ** 2
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
    print("Actual RMSD: " + str(rmsd)[:9] + "     Actual max deviation: " + str(abs_dev)[:9] + "     Worst RMSD: " +
          str(greater_rmsd)[:9] + "     Worst pearson**2: " + str(min(statistic_list_pear))[:9], end="\r")
    return (1 - min(statistic_list_pear)) + (1 - pearson2) + greater_rmsd * 10 + rmsd * 20 + sum(
         statistic_list_pear) / len(statistic_list_pear) + sum(statistic_list_rmsd) / len(statistic_list_rmsd)



if __name__ == "__main__":
    args, logger = arg.settings_argparse()
    warnings.filterwarnings("ignore")
    random.seed(0)
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
        sorted_atomic_types = method.parameters_keys
        for molecule in setm[:number_of_molecules]:
            molecule.symbol_to_number(sorted_atomic_types, method.parameters_type)
        method.make_list_of_lists_of_parameters()
        for molecule in setm[:number_of_molecules]:
            molecule.set_length_correction(method.length_correction)
            try:
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
                                                      args.new_parameters, args.rewriting_with_force, args.chg_output)
        call_stat = "./calculator_charges.py --mode statistics --sdf_input " + args.sdf_input + \
                    "  --charges " + args.right_charges
        call(call_stat, shell=True)
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
        all_atomic_types = []
        for molecule in set_of_molecule:
            all_atomic_types.extend(molecule.symbols(method.parameters_type))
        if args.num_of_parameterized_mol:
            choised_num_of_mol = int(args.num_of_parameterized_mol)
        else:
            choised_num_of_mol = number_of_molecules
        if args.validation:
            high_limit = choised_num_of_mol
            choised_num_of_mol = round(choised_num_of_mol * 0.7)
        set_of_molecule = setm[:choised_num_of_mol]
        atomic_types = para.control_of_missing_atoms(set_of_molecule, method, args.parameters)
        method.set_atomic_types(atomic_types)
        logger.info(colored("Loading molecule data from " + str(args.sdf_input) + " was successful.\n\n\n",
                            "green"))
        method.set_sorted_parameters
        logger.info("Loading charges data from " + str(args.right_charges) + " ...")
        method.load_charges_for_par(args.right_charges)
        method.load_charges_for_par_by_atom_types(args.right_charges, atomic_types)
        logger.info(colored("Loading of charges data was sucessfull. \n\n\n", "green"))
        logger.info(str(number_of_molecules) + " molecules was loaded.\n\n\n")
        input_parameters_list = method.sorted_parameters_values
        if args.method == "EEM":
            bounds = [(0.0001, 2)] * len(input_parameters_list)
        elif args.method == "SAEFM":
            bounds = [(-2, 2)] * len(input_parameters_list)
        elif args.method == "QEq":
            bounds = [(0.0001, 4)] * len(input_parameters_list)
        elif args.method == "SFKEEM":
            bounds = [(0.0001, 4)] * len(input_parameters_list)
        else:
            bounds = [(-1, 1)] * len(input_parameters_list)
        set_of_molecule = setm[:choised_num_of_mol]
        sorted_atomic_types = method.parameters_keys
        for molecule in set_of_molecule:
            molecule.symbol_to_number(sorted_atomic_types, method.parameters_type)
            molecule.set_length_correction(method.length_correction)
        method.make_list_of_lists_of_parameters()
        print("Parameterization running for " + str(choised_num_of_mol) + " molecules ...")
        if args.method_parameterization == "basinhopping":
            res = basinhopping(calculating_charges, input_parameters_list, niter=1, T=0.1, stepsize=0.05)
        elif args.method_parameterization == "diferential_evolution":
            res = differential_evolution(calculating_charges, bounds, maxiter=3)
        else:
            res = minimize(calculating_charges, input_parameters_list, method="SLSQP", bounds=bounds)
            args.method_parameterization = "minimize"
        final_parameters = res.x
        para.writing_new_parameters(args.parameters, args.new_parameters, final_parameters, method)
        logger.info("\n\n\nFile with parameters: \n---------------------------------------")
        logger.info(colored(open(args.new_parameters, 'r').read(), "yellow"))
        logger.info("---------------------------------------")
        if args.validation:
            set_of_molecule = setm[choised_num_of_mol:high_limit]
            for molecule in set_of_molecule:
                molecule.symbol_to_number(sorted_atomic_types, method.parameters_type)
                molecule.set_length_correction(method.length_correction)
        method.make_list_of_lists_of_parameters()
        list_with_data = []
        all_atomic_types_calc = []
        for molecule in set_of_molecule:
            results = method.calculate(molecule)
            list_with_data.append(calc.writing_to_list(molecule, results))
            all_atomic_types_calc.extend(molecule.symbols(method.parameters_type))
        logger.info("Writing calculated charges to " + args.chg_output + " ...")
        calc.writing_to_file(list_with_data, args.chg_output)
        logger.info(colored("Writing to " + args.chg_output + " was successful.\n\n\n", "green"))
        fig = plt.figure(figsize=(11, 9))
        fig_all = fig.add_subplot(111)
        save_fig = args.chg_output[:-4] + ".png"
        args.charges = args.chg_output
        comp.control_if_arguments_files_exist_for_com(args.charges, args.right_charges, save_fig, args.save_fig,
                                                      args.all_mol_to_log, args.rewriting_with_force)
        if args.all_mol_to_log is not None:
            with open(args.all_mol_to_log, "a") as mol_log:
                mol_log.write("name   rmsd   max.dev.   av.dev   pearson**2\n")
        char2 = comp.making_dictionary_with_charges_para(args.charges, all_atomic_types_calc)
        char1 = comp.making_dictionary_with_charges_para(args.right_charges, all_atomic_types)
        list_with_atomic_data, dict_with_atomic_data, dictionary_with_molecular_data, atoms = comp.making_final_list(char1, char2)
        print("Statistics for all atoms:")
        table_for_all_atoms, axis_range, number_of_atoms = comp.statistics_for_all_atoms(list_with_atomic_data,
                                                                                         fitting=True)
        RMSD = table_for_all_atoms[0][0]
        pearson_2 = table_for_all_atoms[0][3]
        print(tabulate(table_for_all_atoms, headers=["RMSD", "max deviation", "average deviation", "pearson**2",
                                                     "num. of atoms"]))
        print("\n\nStatistics for molecules:")
        table_for_all_molecules = comp.statistics_for_molecules(dictionary_with_molecular_data, args.all_mol_to_log)
        print(tabulate(table_for_all_molecules, headers=["RMSD", "max deviation", "average deviation", "pearson**2",
                                                         "num. of molecules"]))
        statistics_data = []
        atoms = sorted(atoms)
        for atom in atoms:
            statistics_data.append(comp.statistics_for_atom_type(atom, dict_with_atomic_data[atom], atoms,
                                                                 fig_all, args.charges, args.right_charges,
                                                                 args.save_fig, axis_range, save_fig[:-4]))
        print("\n\nStatistics for atomic type:")
        if args.method_parameterization:
            method_par = args.method_parameterization
        else:
            method_par = "minimize"
        now = datetime.datetime.now()
        para.write_to_para(args.parameters, args.sdf_input, method_par, choised_num_of_mol, table_for_all_atoms, table_for_all_molecules, statistics_data, now, args.validation)
        print(tabulate(statistics_data, headers=["atomic type", "RMSD", "max deviation", "average deviation",
                                                 "pearson**2", "num. of atoms"]))
        print("\n\n\n")
        if args.make_html:
            make_html.make_html(args.charges[:-4], args.sdf_input, args.method, table_for_all_atoms[0], table_for_all_molecules[0], save_fig, atoms, statistics_data, now.strftime("%Y-%m-%d %H:%M"), args.method_parameterization, choised_num_of_mol, args.validation)
        if args.alarm_after_para:
            alarm.alarm()
        comp.plotting(args.charges, args.right_charges, args.save_fig, save_fig,
                      fig_all, fig, RMSD, pearson_2, axis_range, number_of_atoms)

    elif args.mode == "comparison":
        fig = plt.figure(figsize=(11, 9))
        fig_all = fig.add_subplot(111)
        save_fig = args.charges[:-4] + ".png"
        comp.control_if_arguments_files_exist_for_com(args.charges, args.right_charges, save_fig, args.save_fig,
                                                      args.all_mol_to_log, args.rewriting_with_force)
        if args.all_mol_to_log is not None:
            with open(args.all_mol_to_log, "a") as mol_log:
                mol_log.write("name   rmsd   max.dev.   av.dev   pearson**2\n")
        logger.info("Loading data from " + args.charges + " and " + args.right_charges + " ...")
        char2 = comp.making_dictionary_with_charges(args.charges)
        char1 = comp.making_dictionary_with_charges(args.right_charges)
        list_with_atomic_data, dict_with_atomic_data, dictionary_with_molecular_data, atoms = comp.making_final_list(char1, char2)
        logger.info(
            colored("Loading data from " + args.charges + " and " + args.right_charges + " was sucessfull.\n\n\n",
                    "green"))
        print("Statistics for all atoms:")
        table_for_all_atoms, axis_range, number_of_atoms = comp.statistics_for_all_atoms(list_with_atomic_data,
                                                                                         fitting=True)
        RMSD = table_for_all_atoms[0][0]
        pearson_2 = table_for_all_atoms[0][3]
        print(tabulate(table_for_all_atoms, headers=["RMSD", "max deviation", "average deviation", "pearson**2",
                                                     "num. of atoms"]))
        print("\n\nStatistics for molecules:")
        table_for_all_molecules = comp.statistics_for_molecules(dictionary_with_molecular_data, args.all_mol_to_log)
        print(tabulate(table_for_all_molecules, headers=["RMSD", "max deviation", "average deviation", "pearson**2",
                                                         "num. of molecules"]))
        statistics_data = []
        atoms = sorted(atoms)
        for atom in atoms:
            statistics_data.append(comp.statistics_for_atom_type(atom, dict_with_atomic_data[atom], atoms,
                                                                 fig_all, args.charges, args.right_charges,
                                                                 args.save_fig, axis_range, save_fig[:-4]))
        print("\n\nStatistics for atomic type:")
        print(tabulate(statistics_data, headers=["atomic type", "RMSD", "max deviation", "average deviation",
                                                 "pearson**2", "num. of atoms"]))
        print("\n\n\n")
        comp.plotting(args.charges, args.right_charges, args.save_fig, save_fig,
                      fig_all, fig, RMSD, pearson_2, axis_range, number_of_atoms)

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
                exit("\n")

    elif args.mode == "make_complete_html":
        make_html.make_complete_html(args.verbose)

    elif args.mode == "clusterization":
        clust.control_if_arguments_files_exist_for_clust(args.charges, args.sdf_input)
        logger.info("Checking, if sdf file and chg file are for the same molecules...")
        call_stat = "./calculator_charges.py --mode statistics --sdf_input " + args.sdf_input + \
                    "  --charges " + args.charges
        call(call_stat, shell=True)
        logger.info(colored("Sdf file and chg file are for the same molecules. \n\n\n", "green"))
        logger.info("Loading data... from {} and {}.".format(args.sdf_input, args.charges))
        atom_type = clust.check_type_of_atom(args.atom_type_for_clusterization)
        charges_data, charges_values = clust.loadig_data_from_sdf_and_chg_file(args.sdf_input, args.charges, atom_type,
                                                                          args.atom_type_for_clusterization)
        logger.info(colored("Loading of data was sucessfull. \n\n\n", "green"))
        if args.clusters == "0":
            clust.graph_of_all_charges(args.atom_type_for_clusterization, charges_values, args.fine_of_graph, args.save_fig)
        logger.info("Clusterization of {} clusters...".format(args.clusters))
        counts, average_charges, charges_data = clust.clusterize(charges_values, args.clusters, charges_data)
        logger.info(colored("Clusterization was sucessfull. \n\n\n", "green"))
        clust.print_clusters(args.clusters, counts, average_charges, charges_data)








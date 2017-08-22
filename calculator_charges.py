#!/usr/bin/python

from sys import exit, path, stdout
import argparse
import argcomplete
from numpy import linalg, sqrt
import importlib
from termcolor import colored
from subprocess import call
from scipy.optimize import minimize
from matplotlib import pyplot as plt
from tabulate import tabulate
import warnings
path.insert(0, './modules')
from set_of_molecule import Set_of_molecule
from func_calc import *
from func_para import *
from func_comp import *


def settings_argparse():
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", help="It can be choosen calculation.", required=True, choices=("calculation", "parameterization", "comparison"))
    parser.add_argument("--method", help="Method to calculating charges.")
    parser.add_argument("--sdf_input", help="Sdf file with molecules data.")
    parser.add_argument("--parameters", help="File with parameters.")
    parser.add_argument("--chg_output", help="Output file with charges.")
    parser.add_argument("--right_charges", help="File with charges from for parameterizatin.")
    parser.add_argument("--new_parameters", help="File to save parameters from parameterization.")
    parser.add_argument("--charges1", help="File with charges for comparison.")
    parser.add_argument("--charges2", help="File with charges for comparison")
    parser.add_argument("--save_fig", help="Save figure of comparison.")
    parser.add_argument("--comparison_after_par", action="store_true", help="Comparison after parameterization.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity.")
    argcomplete.autocomplete(parser)
    args = parser.parse_args()
    if args.mode == "calculation":
        if args.method is None or args.sdf_input is None or args.parameters is None or args.chg_output is None:
            parser.error("For calculation must be choosen --method, --sdf_input, --parameters and --chg_output!")
    if args.mode == "parameterization":
        if args.right_charges is None or args.method is None or args.sdf_input is None or args.parameters is None or args.new_parameters is None or args.chg_output is None:
            parser.error("For parameterization must be chossen --right_charges, --method, --sdf_input, --parameters, --chg_output and --new_parameters!")
    if args.mode == "comparison":
        if args.charges1 is None or args.charges2 is None:
            parser.error("For comparison must be chossen --charges1 and --charges2!")
    if args.verbose:
        print("\n\n")
        print("---------------------------------------")
        print("\n\n")
        print(colored("Verbosity turned on.\nMode:" + args.mode + "\n\n\n", "blue"))
    if args.verbose is False:
        print("\n\n\n")


def calculating_charges(list_of_parameters):
    method.load_parameters_from_list(list_of_parameters)
    list_with_results = []
    for molecule in set_of_molecule:
        try:
            results = method.calculate(molecule)
        except linalg.linalg.LinAlgError:
            continue
        list_with_results.extend(results[:-1])
    deviation_list = zip(list_with_results, method.right_charges_for_parametrization)
    rmsd_list = []
    for atom in deviation_list:
        rmsd_list.append((atom[0] - atom[1])**2)
    try:
        rmsd = sqrt((1.0/len(rmsd_list))*sum(rmsd_list))
    except:
        rmsd = 100000.
    stdout.write('Actual RMSD: ' + str(rmsd) + '\r')
    return rmsd


if __name__ == "__main__":
    settings_argparse()
    warnings.filterwarnings("ignore")
    if args.mode == "calculation":
        control_if_arguments_files_exist_for_calc(args.parameters, args.sdf_input, args.chg_output)
        list_with_data = []
        if args.verbose:
            print("Loading molecule data from " + str(args.sdf_input) + ".")
        setm = Set_of_molecule(args.sdf_input)
        number_of_molecules = setm.number_of_molecules
        set_of_molecule = setm.molecules(number_of_molecules)
        if args.verbose:
            print("Loading molecule data from " + str(args.sdf_input) + " was successful.")
            print(str(number_of_molecules) + " molecules was loaded. \n\n\n")
            print("Charges will be calculated by " + str(args.method) + " method.")
        try:
            method = getattr(importlib.import_module("methods"), str(args.method))
        except AttributeError:
            exit(colored("ERROR! Method do not exist or you do not define method!\n\n\n", "red"))
        if args.verbose:
            print("Loading parameters from " + str(args.parameters) + ".")
        method = method()
        method.load_parameters(args.parameters)
        if args.verbose:
            print("File with parameters: \n---------------------------------------")
            with open(args.parameters, 'r') as file:
                print(file.read())
                print("---------------------------------------")
        if method.method_in_parameters != args.method:
            exit(colored("ERROR! These parameters are for method " + str(method.method_in_parameters) + " but you want to calculate charges by method " + str(args.method) + "!\n\n\n", "red"))
        if args.verbose:
            print("Loading of parameters was sucessfull.\n\n\n")
            print("Calculating of charges...\n")
        num_of_err = 0
        for molecule in set_of_molecule:
            try:
                results = method.calculate(molecule)
            except linalg.linalg.LinAlgError:
                print("Molecule " + molecule.name + " is not calculated. Solution do not exist.")
                num_of_err += 1
                continue
            except KeyError:
                print("Molecule " + molecule.name + " contain atoms, for which are not parameters in file: " + args.parameters)
                num_of_err += 1
                continue
            list_with_data.append(writing_to_list(molecule, results))
        if args.verbose:
            print(colored("\n\nCalculation was successful.\n\n\n", "green"))
            print(str(number_of_molecules - num_of_err) + " molecules from " + str(number_of_molecules) + " loaded molecules was calculated. \n\n\n")
        writing_to_file(list_with_data, args.chg_output, args.verbose)

    elif args.mode == "parameterization":
        control_if_arguments_files_exist_for_par(args.right_charges, args.sdf_input, args.parameters, args.new_parameters)
        try:
            method = getattr(importlib.import_module("methods"), str(args.method))
        except AttributeError:
            exit(colored("ERROR! Method do not exist or you do not define method!\n", "red"))
        if args.verbose:
            print("Loading parameters from " + str(args.parameters) + ".")
        method = method()
        method.load_parameters(args.parameters)
        if args.verbose:
            print("File with parameters: \n---------------------------------------")
            with open(args.parameters, 'r') as file:
                print(file.read())
                print("---------------------------------------")
        if method.method_in_parameters != args.method:
            exit(colored("ERROR! These parameters are for method " + str(method.method_in_parameters) + " but you want to calculate charges by method " + str(args.method) + "!\n\n\n", "red"))
        if args.verbose:
            print("Loading of parameters was sucessfull.\n\n\n")
        if args.verbose:
            print("Loading molecule data from " + str(args.sdf_input) + ".")
        setm = Set_of_molecule(args.sdf_input)
        number_of_molecules = setm.number_of_molecules
        set_of_molecule = setm.molecules(number_of_molecules)
        control_of_missing_atoms(set_of_molecule, method, args.parameters)
        if args.verbose:
            print("Loading molecule data from " + str(args.sdf_input) + " was successful.")
            print(str(number_of_molecules) + " molecules was loaded. \n\n\n")
        method.set_sorted_parameters()
        method.load_charges_for_par(args.right_charges)
        input_parameters_list = method.parameters
        for x in [1,25,100]:
            if x >= number_of_molecules:
                break
            print("Parameterization running for " + str(x) + " molecule.\n\n\n")
            set_of_molecule = setm.molecules(x)
            res = minimize(calculating_charges, input_parameters_list, method="SLSQP", options={'disp': True})
            print("\n\n\n")
            input_parameters_list = res.x



        set_of_molecule = setm.molecules(number_of_molecules)

        print("Parameterization running for " + str(number_of_molecules) + " molecule.\n\n\n")
        res = minimize(calculating_charges, input_parameters_list, method="SLSQP", options={'disp': True})




        writing_new_parameters(args.parameters, args.new_parameters, res, method)
        if args.verbose:
            print("\n\nFile with new parameters: \n---------------------------------------")
            with open(args.new_parameters, 'r') as file:
                print(file.read())
                print("---------------------------------------\n\n")
                print(colored("Parameterization was successful.\n\n\n", "green"))
        if args.comparison_after_par:
            call1 = "./calculator_charges.py --mode calculation --parameters " + args.new_parameters + " --sdf_input " + args.sdf_input + " --method " + args.method + " --chg_output " + args.chg_output
            call2 = "./calculator_charges.py --mode comparison --charges1 " + args.chg_output + " --charges2 " + args.right_charges
            if args.verbose:
                call1 = call1 + " -v"
                call2 = call2 + " -v"
            if args.save_fig:
                call2 = call2 + " --save_fig " + args.save_fig
            exit_value = call(call1, shell=True)
            if exit_value == 0:
                call(call2, shell=True)

    elif args.mode == "comparison":
        control_if_arguments_files_exist_for_com(args.charges1, args.charges2, args.save_fig)
        char1 = making_dictionary_with_charges(args.charges1)
        char2 = making_dictionary_with_charges(args.charges2)
        list_with_atomic_data, dictionary_with_molecular_data, atoms = making_final_list(char1, char2)
        print("Statistics for all atoms")
        table_for_all_atoms, number_of_atoms = statistics_for_all_atoms(list_with_atomic_data)
        print(tabulate(table_for_all_atoms, headers=["RMSD", "max deviation", "average deviation", "pearson"]))
        print("Number of all atoms is: " + str(number_of_atoms))
        print("")
        print("")
        print("Statistics for molecules:")
        table_for_all_molecules, number_of_molecules = statistics_for_molecules(dictionary_with_molecular_data)
        print(tabulate(table_for_all_molecules, headers=["RMSD", "max deviation", "average deviation", "pearson"]))
        print("Number of all molecules is: " + str(number_of_molecules))
        statistics_data = []
        fig = plt.figure(figsize=(14, 12))
        for atom in atoms:
            statistics_data.append(statistics_for_atom_type(atom, list_with_atomic_data, atoms))
        print("")
        print("")
        print("Statistics for atomic type:")
        print(
            tabulate(statistics_data, headers=["atomic type", "RMSD", "max deviation", "average deviation", "pearson"]))
        print("")
        print("")
        plotting(args.charges1, args.charges2, args.save_fig)

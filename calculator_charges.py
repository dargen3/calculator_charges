#!/usr/bin/python
from sys import exit, stdout, stdin
import argparse
import argcomplete
from numpy import linalg, sqrt
import importlib
from termcolor import colored
from subprocess import call, check_call
from scipy.optimize import minimize, basinhopping, differential_evolution
from matplotlib import pyplot as plt
from tabulate import tabulate
import warnings
from modules.set_of_molecule import Set_of_molecule
import modules.func_calc as calc
import modules.func_para as para
import modules.func_comp as comp
import modules.func_stat as stat
from numba import jit
from time import time



def settings_argparse():
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", help="It can be choosen calculation, parameterization, comparison or statistics.", required=True, choices=("calculation", "parameterization", "comparison", "statistics"))
    parser.add_argument("--method", help="Method to calculating charges.")
    parser.add_argument("--sdf_input", help="Sdf file with molecules data.")
    parser.add_argument("--parameters", help="File with parameters.")
    parser.add_argument("--chg_output", help="Output file with charges.")
    parser.add_argument("--right_charges", help="File with charges from for parameterizatin.")
    parser.add_argument("--new_parameters", help="File to save parameters from parameterization.")
    parser.add_argument("--charges", help="File with charges for comparison.")
    parser.add_argument("--num_of_parameterized_mol", help="Only first N molecule will be parameterized.")
    parser.add_argument("--save_fig", help="Save figure of comparison.")
    parser.add_argument("--save_setm", help="Not for using. Its for intern usage of script.", action="store_true")
    parser.add_argument("--atom_types_for_para", help="Only for SAEFM!", nargs="+")
    parser.add_argument("--method_parameterization", help="It can be choosen minimize, basinhopping or diferential_evolution.", choices=("minimize", "basinhopping", "diferential_evolution"))
    parser.add_argument("--comparison_after_par", action="store_true", help="Comparison after parameterization.")
    parser.add_argument("--alarm_after_para", action="store_true", help="Alarm after parameterization. You need instalated rhythmbox!")
    parser.add_argument("--statistic_on_end_of_para", help="Statistics data from parameterization is added to end of parameters file.")
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
        if args.charges is None or args.right_charges is None:
            parser.error("For comparison must be chossen --charges and --right_charges!")
    if args.mode == "statistics":
        if args.charges is None or args.sdf_input is None:
            parser.error("For statistics must be chossen --charges and --sdf_input")
    if args.verbose:
        print("\n\n")
        print("---------------------------------------")
        print("\n\n")
        print(colored("Verbosity turned on.\nMode:" + args.mode + "\n\n\n", "blue"))
    if args.verbose is False:
        print("\n\n\n")


@jit(nopython=True)
def statistics(deviation_list):
    rmsd_list = []
    abs_dev = 0
    for atom in deviation_list:
        rmsd_list.append((atom[0] - atom[1])**2)
        dev = abs(atom[0] - atom[1])
        if dev > abs_dev:
            abs_dev = dev
    return rmsd_list, abs_dev


def calculating_charges(list_of_parameters):
    method.load_parameters_from_list(list_of_parameters)
    list_with_results = []
    for molecule in set_of_molecule:
        try:
            if args.method[0:2] == "SA":
                results = method.calculate(molecule, method._get_parameterized_atom_type)
            else:
                results = method.calculate(molecule)
        except linalg.linalg.LinAlgError:
            continue
        list_with_results.extend(results)
    if method._get_parameterized_atom_type is False:
        deviation_list = zip(list_with_results, method._right_charges_for_parametrization)
    else:
        deviation_list = zip(list_with_results, method.right_charges_for_parameterization_by_atom_types(method._get_parameterized_atom_type))
    rmsd_list, abs_dev = statistics(deviation_list)
    rmsd = sqrt((1.0 / len(rmsd_list)) * sum(rmsd_list))
    stdout.write("Actual RMSD: " + str(rmsd) + "     Actual max deviation: " + str(abs_dev) + "\r")
    return rmsd


if __name__ == "__main__":
    settings_argparse()
    warnings.filterwarnings("ignore")
    if args.mode == "calculation":
        calc.control_if_arguments_files_exist_for_calc(args.parameters, args.sdf_input, args.chg_output)
        list_with_data = []
        if args.verbose:
            print("Loading molecule data from " + str(args.sdf_input) + ".")
        setm = Set_of_molecule(args.sdf_input)
        number_of_molecules = setm._number_of_molecules
        set_of_molecule = setm.molecules(number_of_molecules)
        if args.verbose:
            print("Loading molecule data from " + str(args.sdf_input) + " was successful.")
            print(str(number_of_molecules) + " molecules was loaded. \n\n\n")
            print("Charges will be calculated by " + str(args.method) + " method.")
        try:
            method = getattr(importlib.import_module("modules.methods"), str(args.method))
        except AttributeError:
            exit(colored("ERROR! Method do not exist or you do not define method!\n\n\n", "red"))
        if args.verbose:
            print("Loading parameters from " + str(args.parameters) + ".")
        method = method()
        method.load_parameters(args.parameters)
        if args.verbose:
            print("File with parameters: \n---------------------------------------")
            with open(args.parameters, 'r') as file:
                print(colored(file.read(), "yellow"))
                print("---------------------------------------")
        if method._method_in_parameters != args.method:
            exit(colored("ERROR! These parameters are for method " + str(method._method_in_parameters) + " but you want to calculate charges by method " + str(args.method) + "!\n\n\n", "red"))
        if args.verbose:
            print("Loading of parameters was sucessfull.\n\n\n")
            print("Calculating of charges...\n")
        method.set_parameters_type()
        method.set_getting_parameters("calculation")
        num_of_err = 0
        if args.method[0:2] == "SA":
            atomic_types = para.control_of_missing_atoms(set_of_molecule, method, args.parameters)
        for molecule in set_of_molecule:
            molecule.set_length_correction(method._length_correction)
            try:
                if args.method[0:2] == "SA":
                    results = method.calculate(molecule, atomic_types)
                else:
                    results = method.calculate(molecule)
            except linalg.linalg.LinAlgError:
                print("Molecule " + molecule._name + " is not calculated. Solution do not exist.")
                num_of_err += 1
                continue
            except KeyError:
                print("Molecule " + molecule._name + " contain atoms, for which are not parameters in file: " + args.parameters)
                num_of_err += 1
                continue
            list_with_data.append(calc.writing_to_list(molecule, results))
        if args.verbose:
            print(colored("\n\nCalculation was successful.\n\n\n", "green"))
            print(str(number_of_molecules - num_of_err) + " molecules from " + str(number_of_molecules) + " loaded molecules was calculated. \n\n\n")
        calc.writing_to_file(list_with_data, args.chg_output, args.verbose)
    elif args.mode == "parameterization":
        start_time = time()
        para.control_if_arguments_files_exist_for_par(args.right_charges, args.sdf_input, args.parameters, args.new_parameters)
        try:
            check_call("./calculator_charges.py --mode statistics -v --save_setm --sdf_input " + args.sdf_input + "  --charges " + args.right_charges, shell=True)
        except:
            print("asdfasdf")
            exit()
        print("---------------------------------------\n\n\n")
        try:
            method = getattr(importlib.import_module("modules.methods"), str(args.method))
        except AttributeError:
            exit(colored("ERROR! Method do not exist or you do not define method!\n", "red"))
        if args.verbose:
            print("Loading parameters from " + str(args.parameters) + ".")
        method = method()
        method.load_parameters(args.parameters)
        if args.verbose:
            print("File with parameters: \n---------------------------------------")
            with open(args.parameters, 'r') as file:
                print(colored(file.read(), "yellow"))
                print("---------------------------------------")
        if method._method_in_parameters != args.method:
            exit(colored("ERROR! These parameters are for method " + str(method._method_in_parameters) + " but you want to calculate charges by method " + str(args.method) + "!\n\n\n", "red"))
        method.set_parameters_type()
        method.set_getting_parameters("parameterization")
        if args.verbose:
            print("Loading of parameters was sucessfull.\n\n\n")
            print("Loading molecule data from " + str(args.sdf_input) + ".")
        setm = Set_of_molecule(args.sdf_input)
        number_of_molecules = setm._number_of_molecules
        set_of_molecule = setm.molecules(number_of_molecules)
        for molecule in set_of_molecule:
            molecule.set_length_correction(method._length_correction)
        if args.num_of_parameterized_mol:
            chooised_num_of_mol = int(args.num_of_parameterized_mol)
        else:
            chooised_num_of_mol = number_of_molecules
        if args.method[0:2] == "SA":
            set_of_molecule = setm.molecules(chooised_num_of_mol)
            a_types = para.control_of_missing_atoms(set_of_molecule, method, args.parameters)
            if args.atom_types_for_para:
                for x in args.atom_types_for_para:
                    if x not in a_types:
                        exit(colored("ERROR! In chooised molecules are no atoms with symbol " + x,"red"))
                atomic_types = args.atom_types_for_para
            else:
                atomic_types = a_types
            if args.verbose:
                print("Loading molecule data from " + str(args.sdf_input) + " was successful.\n\n\n")
            global_sorted_parameters = method.set_sorted_parameters()
            if args.verbose:
                print("Loading charges data from " + str(args.right_charges) + ".")
            num_of_charges = method.load_charges_for_par_by_atom_types(args.right_charges, atomic_types)
            if args.verbose:
                print("Loading of charges data was sucessfull. \n\n\n")
            if args.verbose:
                print(str(number_of_molecules) + " molecules was loaded.\n\n\n")
            global_input_parameters_list = method.parameters
            final_parameters = global_input_parameters_list
            method.load_charges_for_par_with_molecules(args.right_charges)
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
                print("Parameterization running for atom: " + atomic_type)
                bounds = [(-1, 5)] * len(input_parameters_list)
                if args.method_parameterization == "basinhopping":
                    res = basinhopping(calculating_charges, input_parameters_list, niter=1, T=0.1, stepsize=0.05)
                elif args.method_parameterization == "diferential_evolution":
                    res = differential_evolution(calculating_charges, bounds, maxiter=1)
                else:
                    res = minimize(calculating_charges, input_parameters_list, method="SLSQP", options={'disp': True}, bounds=bounds)
                for x in range(len(res.x)):
                    final_parameters[global_sorted_parameters.index(sorted_parameters[x])] = res.x[x]
                print("\n\n")
            method.set_global_sorted_parameters_keys()
        else:
            set_of_molecule = setm.molecules(chooised_num_of_mol)
            para.control_of_missing_atoms(set_of_molecule, method, args.parameters)
            if args.verbose:
                print("Loading molecule data from " + str(args.sdf_input) + " was successful.\n\n\n")
            method.set_sorted_parameters()
            if args.verbose:
                print("Loading charges data from " + str(args.right_charges) + ".")
            method.load_charges_for_par(args.right_charges)
            if args.verbose:
                print("Loading of charges data was sucessfull. \n\n\n")
            if args.verbose:
                print(str(number_of_molecules) + " molecules was loaded atoms.\n\n\n")
            input_parameters_list = method.parameters
            bounds = [(-2, 2)] * len(input_parameters_list)
            for x in [10,100,1000]:
                if x >= chooised_num_of_mol:
                    break
                print("Parameterization running for " + str(x) + " molecules.\n\n\n")
                set_of_molecule = setm.molecules(x)
                if args.method_parameterization == "basinhopping":
                    res = basinhopping(calculating_charges, input_parameters_list, niter=1, T=0.1, stepsize=0.05)
                elif args.method_parameterization == "diferential_evolution":
                    res = differential_evolution(calculating_charges, bounds, maxiter=1)
                else:
                    res = minimize(calculating_charges, input_parameters_list, method="SLSQP", options={'disp': True}, bounds=bounds)
                print("\n\n\n")
                input_parameters_list = res.x
            set_of_molecule = setm.molecules(chooised_num_of_mol)
            print("Parameterization running for " + str(number_of_molecules) + " molecules.\n\n\n")
            if args.method_parameterization == "basinhopping":
                res = basinhopping(calculating_charges, input_parameters_list, niter=1, T=0.1, stepsize=0.05)
            elif args.method_parameterization == "diferential_evolution":
                res = differential_evolution(calculating_charges, bounds, maxiter=1)
            else:
                res = minimize(calculating_charges, input_parameters_list, method="SLSQP", options={'disp': True}, bounds=bounds)
            final_parameters = res.x
        para.writing_new_parameters(args.parameters, args.new_parameters, final_parameters, method)
        if args.verbose:
            print("\n\nFile with new parameters: \n---------------------------------------")
            with open(args.new_parameters, 'r') as file:
                print(colored(file.read(), "yellow"))
                print("---------------------------------------\n\n")
                print(colored("Parameterization was successful.\n\n\n", "green"))
        if args.alarm_after_para:
            try:
                check_call("play alarm/alarm.mp3 2> /dev/null", shell=True)
            except KeyboardInterrupt:
                pass
            except:
                print(colored("You cannot run mp3 by bash command play!", "red"))
        if args.comparison_after_par:
            call1 = "./calculator_charges.py --mode calculation --parameters " + args.new_parameters + " --sdf_input " + args.sdf_input + " --method " + args.method + " --chg_output " + args.chg_output
            call2 = "./calculator_charges.py --mode comparison --charges " + args.chg_output + " --right_charges " + args.right_charges + " --statistic_on_end_of_para " + args.new_parameters + " --sdf_input " + args.sdf_input
            if args.verbose:
                call1 = call1 + " -v"
                call2 = call2 + " -v"
            if args.save_fig:
                call2 = call2 + " --save_fig " + args.save_fig
            exit_value = call(call1, shell=True)
            if exit_value == 0:
                call(call2, shell=True)
            end_time = time()
            time_taken = end_time - start_time
            hours, rest = divmod(time_taken, 3600)
            minutes, seconds = divmod(rest, 60)
            seconds = round(seconds)
            print("\n\nParameterization take " + str(hours)[:-2] + " hours " + str(minutes)[:-2] + " minutes " + str(seconds)[:-2] + " seconds.\n\n\n")

    elif args.mode == "comparison":
        comp.control_if_arguments_files_exist_for_com(args.charges, args.right_charges, args.save_fig)
        char2 = comp.making_dictionary_with_charges(args.charges)
        char1 = comp.making_dictionary_with_charges(args.right_charges)
        list_with_atomic_data, dictionary_with_molecular_data, atoms = comp.making_final_list(char1, char2)
        print("Statistics for all atoms:")
        table_for_all_atoms, number_of_atoms = comp.statistics_for_all_atoms(list_with_atomic_data)
        print(tabulate(table_for_all_atoms, headers=["RMSD", "max deviation", "average deviation", "pearson", "num. of atoms"]))
        if args.statistic_on_end_of_para:
            with open(args.statistic_on_end_of_para, "a") as parameters:
                parameters.write("\n\n\nParameterized set of molecules: " + args.sdf_input + "\n")
                if args.method_parameterization:
                    method_par = args.method_parameterization
                else:
                    method_par = "minimize"
                parameters.write("Method of parameterization: " + method_par + "\n\n")
                parameters.write("Statistics for all atoms:\n")
                parameters.write(tabulate(table_for_all_atoms, headers=["RMSD", "max deviation", "average deviation", "pearson", "num. of atoms"]))
                parameters.write("\n\n")
        print("")
        print("")
        print("Statistics for molecules:")
        table_for_all_molecules = comp.statistics_for_molecules(dictionary_with_molecular_data)
        print(tabulate(table_for_all_molecules, headers=["RMSD", "max deviation", "average deviation", "pearson", "num. of molecules"]))
        if args.statistic_on_end_of_para:
            with open(args.statistic_on_end_of_para, "a") as parameters:
                parameters.write("Statistics for molecules:\n")
                parameters.write(tabulate(table_for_all_molecules, headers=["RMSD", "max deviation", "average deviation", "pearson", "num. of molecules"]))
                parameters.write("\n\n")
        statistics_data = []
        fig = plt.figure(figsize=(14, 12))
        for atom in atoms:
            statistics_data.append(comp.statistics_for_atom_type(atom, list_with_atomic_data, atoms))
        print("")
        print("")
        print("Statistics for atomic type:")
        print(tabulate(statistics_data, headers=["atomic type", "RMSD", "max deviation", "average deviation", "pearson", "num. of atoms"]))
        if args.statistic_on_end_of_para:
            with open(args.statistic_on_end_of_para, "a") as parameters:
                parameters.write("Statistics for atomic type:\n")
                parameters.write(tabulate(statistics_data, headers=["atomic type", "RMSD", "max deviation", "average deviation", "pearson", "num. of atoms"]))
                parameters.write("\n\n")
        print("")
        print("")
        comp.plotting(args.charges, args.right_charges, args.save_fig)
    elif args.mode == "statistics":
        charges_data, sdf_data, charges_names, sdf_names, setm, number_of_lines = stat.is_the_same(args.charges, args.sdf_input)
        if args.save_setm:
            with open("setm_file.helped.txt", "w") as setm_file:
                setm_file.write(str(setm))
        if charges_data == sdf_data and sdf_names == charges_names:
                print(colored("Sdf file is for the same molecules like charge file!\n\n\n", "green"))
        else:
            if charges_data == sdf_data and sdf_names != charges_names:
                print(colored("Sdf file is for the same molecules like charge file, but molecules is in different order!\n\n\n", "red"))
            else:
                print(colored("Sdf file is for diferent set of molecules like charges file!", "red"))
            print("If you want to create new sdf and chg file only with together molecules write yes. Else press enter.")
            decision = stdin.readline().rstrip('\n')
            if decision == "yes":
                stat.create_new_files(sdf_names, charges_names, args.sdf_input, args.charges, setm, number_of_lines)
                exit(colored("Cresting new files was successful!", "green"))
            else:
                exit(colored("\n\n Statistics was successful!", "green"))




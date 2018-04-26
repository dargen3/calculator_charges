#!/usr/bin/env python3
import warnings
from numpy import random
from modules.calculation import calculate_charges
from modules.parameterization import parameterize
from modules.comparison import comparison
from modules.statistics import statistics
from modules.clusterization import clusterize
from modules.make_html import make_complete_html
from modules.set_argparse import settings_argparse
from modules.parameterization_find_args import find_argumets_and_parameterize
from modules.set_of_molecule import Set_of_molecule

if __name__ == "__main__":
    args, logger = settings_argparse()
    warnings.filterwarnings("ignore")
    # random.seed(0)
    if args.mode == "calculation":
        calculate_charges(args.parameters, args.sdf_input, args.chg_output, args.rewriting_with_force, args.method,
                          logger)

    elif args.mode == "parameterization":
        parameterize(args.method, args.parameters, args.sdf_input, args.num_of_parameterized_mol, args.validation,
                     args.right_charges, args.method_parameterization, args.new_parameters, args.chg_output,
                     args.all_mol_to_log, logger, args.rewriting_with_force, args.save_fig, args.make_html,
                     args.cpu)

    elif args.mode == "comparison":
        comparison(args.sdf_input, args.charges, args.right_charges, args.save_fig, args.all_mol_to_log, args.rewriting_with_force,
                   logger)

    elif args.mode == "statistics":
        statistics(args.charges, args.sdf_input, logger)

    elif args.mode == "make_complete_html":
        make_complete_html(args.verbose)

    elif args.mode == "clusterization":
        clusterize(args.charges, args.sdf_input, logger, args.atom_type_for_clusterization, args.clusters,
                   args.save_fig, args.fine_of_graph)

    elif args.mode == "parameterization_find_args":
        find_argumets_and_parameterize(args.path, logger, args.mode, args.method_parameterization, args.num_of_parameterized_mol, args.cpu, args.time)

    elif args.mode == "parameterization_send_meta":
        find_argumets_and_parameterize(args.path, logger, args.mode, args.method_parameterization, args.num_of_parameterized_mol, args.cpu, args.time)

    elif args.mode == "set_of_molecule_info":
        Set_of_molecule(args.sdf_input).statistics_data()

    elif args.mode == "atoms_data":
        with open(args.right_charges, "r") as right_charges_file:
            list_with_right_charges = []
            for line in right_charges_file:
                ls = line.split()
                if len(ls) not in {0, 1, 3}:
                    exit(colored("File " + file + " with right charges is wrong! \n", "red"))
                if len(ls) == 3:
                    list_with_right_charges.append(float(ls[2]))
        setm = Set_of_molecule(args.sdf_input)
        set_of_molecules = setm[:len(setm)]
        count = 0
        for molecule in set_of_molecules:
            for atom in molecule._atoms:
                atom.charge = list_with_right_charges[count]
                count += 1


        electronegativity = {"C":2.55, "H":2.2, "O":3.44, "N":3.04, "S":2.58}
        for molecule in set_of_molecules:
            for index, atom in enumerate(molecule._atoms):
                atom.highest_bond = molecule._symbol_gravity[index][-1]
                atom.atoms_data = {2:[0,0,0],3:[0,0,0],4:[0,0,0],5:[0,0,0],6:[0,0,0],7:[0,0,0],8:[0,0,0],9:[0,0,0],10:[0,0,0],11:[0,0,0],12:[0,0,0]}
                for i,a in enumerate(molecule._atoms):
                    if i == index:
                        continue
                    for key in atom.atoms_data:
                        if molecule.matrix_of_distance[i][index] < key:
                            atom.atoms_data[key][0] += (electronegativity[a.symbol] - electronegativity[atom.symbol])/molecule.matrix_of_distance[index][i]
                            atom.atoms_data[key][1] += 1
                            atom.atoms_data[key][2] += (a.charge - atom.charge)/molecule.matrix_of_distance[index][i]
        print("type,charge,highest_bond,el2,chg2,el3,chg3,el4,chg4,el5,chg5,el6,chg6,el7,chg7,el8,chg8,el9,chg9,el10,chg10,el11,chg11,el12,chg12")
        for molecule in set_of_molecules:
            for a in molecule._atoms:
                print(a.symbol, end=",")
                print(a.charge, end=",")
                print(a.highest_bond, end=",")
                for dist in range(2,13):
                    print(a.atoms_data[dist][0]/a.atoms_data[dist][1], end=",")
                    if dist != 12:
                        print(a.atoms_data[dist][2]/a.atoms_data[dist][1], end=",")
                    else:
                        print(a.atoms_data[dist][2]/a.atoms_data[dist][1])







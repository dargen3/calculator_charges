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
        electronegativity = {"C":2.55, "H":2.2, "O":3.44, "N":3.04, "S":2.58}
        print("type,charge,highest_bond,el_BA,chg_BA,el_BBA,chg_BBA,total_chg,el_dist,el_dist_cor,total_chg_cor")
        for molecule in set_of_molecules:
            for index, atom in enumerate(molecule._atoms):
                atom.charge = list_with_right_charges[count]
                count += 1
                atom.highest_bond = molecule._symbol_gravity[index][-1]
                atom.num_of_BA = len(molecule._bonded_atoms[index+1])
                atom.BA = molecule._bonded_atoms[index+1]
                atom.num_of_BBA = len(molecule._bonded_bonded_atoms[index+1])
                atom.BBA = molecule._bonded_bonded_atoms[index+1]
                atom.el_BA = 0
                atom.el_BBA = 0
                atom.el_dist = 0
                atom.el_dist_cor = 0
                for aaatom in molecule._bonded_atoms[index+1]:
                    atom.el_BA += electronegativity[molecule._atoms[aaatom-1].symbol]
                    atom.el_dist += electronegativity[molecule._atoms[aaatom-1].symbol]/molecule.matrix_of_distance[index][aaatom-1]
                    atom.el_dist_cor += (electronegativity[molecule._atoms[aaatom-1].symbol] - electronegativity[atom.symbol])/molecule.matrix_of_distance[index][aaatom-1]
                for aaatom in molecule._bonded_bonded_atoms[index + 1]:
                    atom.el_BBA += electronegativity[molecule._atoms[aaatom - 1].symbol]
                    atom.el_dist += electronegativity[molecule._atoms[aaatom-1].symbol]/molecule.matrix_of_distance[index][aaatom-1]
                    atom.el_dist_cor += (electronegativity[molecule._atoms[aaatom - 1].symbol] - electronegativity[atom.symbol]) / molecule.matrix_of_distance[index][aaatom - 1]

        for molecule in set_of_molecules:
            for index, atom in enumerate(molecule._atoms):
                atom.chg_BA = 0
                atom.chg_BBA = 0
                atom.total_chg = 0
                atom.total_chg_cor = 0
                for aaaa in atom.BA:
                    atom.chg_BA += molecule._atoms[aaaa-1].charge
                    atom.total_chg += molecule._atoms[aaaa-1].charge/molecule.matrix_of_distance[index][aaaa-1]
                    atom.total_chg_cor += (molecule._atoms[aaaa-1].charge - atom.charge)/molecule.matrix_of_distance[index][aaaa-1]
                for bbbb in atom.BBA:
                    atom.chg_BBA += molecule._atoms[bbbb-1].charge
                    atom.total_chg_cor += (molecule._atoms[aaaa-1].charge - atom.charge)/molecule.matrix_of_distance[index][aaaa-1]
                    atom.total_chg += molecule._atoms[aaaa-1].charge/molecule.matrix_of_distance[index][bbbb-1]
                a = atom
                print("{},{},{},{},{},{},{},{},{},{},{}".format(a.symbol, a.charge, a.highest_bond, round(a.el_BA/a.num_of_BA,4), round(a.chg_BA/a.num_of_BA, 4), round(a.el_BBA/a.num_of_BBA,4), round(a.chg_BBA/a.num_of_BBA, 4), round(a.total_chg/(a.num_of_BBA+a.num_of_BA), 4), round(a.el_dist/(a.num_of_BBA+a.num_of_BA), 4), round((a.el_dist_cor/(a.num_of_BBA+a.num_of_BA)-electronegativity[a.symbol]), 4), round((a.total_chg_cor/(a.num_of_BBA+a.num_of_BA))-a.charge, 4)))
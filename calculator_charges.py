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
    random.seed(0)
    if args.mode == "calculation":
        calculate_charges(args.parameters, args.sdf_input, args.chg_output, args.rewriting_with_force, args.method,
                          logger)

    elif args.mode == "parameterization":
        parameterize(args.method, args.parameters, args.sdf_input, args.num_of_parameterized_mol, args.validation,
                     args.right_charges, args.method_parameterization, args.new_parameters, args.chg_output,
                     args.all_mol_to_log, logger, args.rewriting_with_force, args.save_fig, args.make_html,
                     args.alarm_after_para)

    elif args.mode == "comparison":
        comparison(args.charges, args.right_charges, args.save_fig, args.all_mol_to_log, args.rewriting_with_force,
                   logger)

    elif args.mode == "statistics":
        statistics(args.charges, args.sdf_input, logger)

    elif args.mode == "make_complete_html":
        make_complete_html(args.verbose)

    elif args.mode == "clusterization":
        clusterize(args.charges, args.sdf_input, logger, args.atom_type_for_clusterization, args.clusters,
                   args.save_fig, args.fine_of_graph)

    elif args.mode == "parameterization_find_args":
        find_argumets_and_parameterize(args.path, logger, args.mode, args.method_parameterization, args.num_of_parameterized_mol)

    elif args.mode == "parameterization_send_meta":
        find_argumets_and_parameterize(args.path, logger, args.mode, args.method_parameterization, args.num_of_parameterized_mol)

    elif args.mode == "set_of_molecule_info":
        Set_of_molecule(args.sdf_input).statistics_data()

    elif args.mode == "experimental":
        set_of_molecule = Set_of_molecule(args.sdf_input)
        atom_types = []
        for molecule in set_of_molecule[:len(set_of_molecule)]:
            atom_types.extend(molecule.get_all_atom_types())
        from pprint import pprint
        pprint(sorted(set(atom_types)))
        for x in sorted(set(atom_types)):
            print("{} 1".format(x))
    print("\n")








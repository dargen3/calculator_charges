import argparse
import logging
import argcomplete
from termcolor import colored


def settings_argparse():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", help="It can be chosen calculation, parameterization, comparison, statistics, "
                                       "make_complete_html, parrameterization_find_args, parameterization_meta"
                                       "set_of_molecule_info and clusterization.",
                        required=True, choices=("calculation", "parameterization", "comparison", "statistics",
                                                "make_complete_html", "clusterization", "parameterization_find_args",
                                                "set_of_molecule_info", "parameterization_send_meta", "atoms_data"))
    parser.add_argument("--method", help="Method to calculating charges or for parameterization.")
    parser.add_argument("--path", help="For mode --parameterization_find_args")
    parser.add_argument("--sdf_input", help="Sdf file with molecules data.")
    parser.add_argument("--parameters", help="File with parameters.")
    parser.add_argument("--time", help="How time take calculation in META (hours).", default=100)
    parser.add_argument("--fine_of_graph", help="For clusterization --cluster 0. Default is 0.1", type=float)
    parser.add_argument("--all_mol_to_log", help="For comparison. Results of all molecules are saved into log file.")
    parser.add_argument("--chg_output", help="Output chg file with charges.")
    parser.add_argument("--right_charges", help="File with charges from QM.")
    parser.add_argument("--new_parameters", help="File to save parameters from parameterization.")
    parser.add_argument("--charges", help="File with charges.")
    parser.add_argument("--cpu", help="Number of cpu for parameterization", type=int, default=1)
    parser.add_argument("--clusters", help="Number of cluster for mode clusterization.", type=int)
    parser.add_argument("--atom_type_for_clusterization", help="Atom type for clusterizing in format "
                                                               "atom~highest_bond (C~2).")
    parser.add_argument("--num_of_parameterized_mol", help="Only first N molecule will be parameterized.", type=int)
    parser.add_argument("--save_fig", action="store_true", help="Save figures of comparison.")
    parser.add_argument("--validation", action="store_true",
                        help="From set of molecules will be 70% used for parameterization and 30% for validation.")
    parser.add_argument("--method_parameterization",
                        help="It can be chosen minimize or guided_minimization.",
                        choices=("minimize", "guided_minimization"), default="guided_minimization")
    parser.add_argument("--make_html", action="store_true", help="Make html after comparison. Only for "
                                                                 "parameterization. You must choise --save_fig too.")
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
            parser.error("For statistics must be choisen --charges and --sdf_input!")
    elif args.mode == "clusterization":
        if args.charges is None or args.sdf_input is None or args.clusters is None:
            parser.error("For clusterization must be choisen --charges, --sdf_input and --clusters!")
    elif args.mode == "statistics":
        if args.path is None:
            parser.error("For para_find_args must be choisen --path!")
    elif args.mode == "set_of_molecule_info":
        if args.sdf_input is None:
            parser.error("For set of molecule info must be choisen --sdf_input!")
    elif args.mode == "parameterization_find_args":
        if args.path is None:
            parser.error("For parameterization_find_args must be choisen --path!")
    elif args.mode == "parameterization_send_meta":
        if args.path is None:
            parser.error("For parameterization_send_meta must be choisen --path and --method_parameterization!")
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
    logger.info(colored("Verbosity turned on.\nMode: " + args.mode + "\n\n\n", "blue"))
    print("\n")
    return args, logger

import argparse
import logging
import argcomplete
from termcolor import colored


def settings_argparse():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", help="It can be chosen calculation, parameterization, comparison, statistics, make_complete_html and clusterization.",
                        required=True, choices=("calculation", "parameterization", "comparison", "statistics", "make_complete_html", "clusterization"))
    parser.add_argument("--method", help="Method to calculating charges or for parameterization.")
    parser.add_argument("--sdf_input", help="Sdf file with molecules data.")
    parser.add_argument("--parameters", help="File with parameters.")
    parser.add_argument("--fine_of_graph", help="For clusterization --cluster 0. Default is 0.1")
    parser.add_argument("--all_mol_to_log", help="For comparison. Results of all molecules are saved into log file.")
    parser.add_argument("--chg_output", help="Output chg file with charges.")
    parser.add_argument("--right_charges", help="File with charges from QM.")
    parser.add_argument("--new_parameters", help="File to save parameters from parameterization.")
    parser.add_argument("--charges", help="File with charges.")
    parser.add_argument("--clusters", help="Number of cluster for mode clusterization.")
    parser.add_argument("--atom_type_for_clusterization", help="Atom type for clusterizing in format atom~highest_bond (C~2).")
    parser.add_argument("--num_of_parameterized_mol", help="Only first N molecule will be parameterized.")
    parser.add_argument("--save_fig", action="store_true", help="Save figures of comparison.")
    parser.add_argument("--validation", action="store_true",
                        help="From set of molecules will be 70% used for parameterization and 30% for validation.")
    parser.add_argument("--method_parameterization",
                        help="It can be chosen minimize, basinhopping or diferential_evolution.",
                        choices=("minimize", "basinhopping", "diferential_evolution"))
    parser.add_argument("--alarm_after_para", action="store_true",
                        help="Alarm after parameterization. You need instalated bash's play!")
    parser.add_argument("--make_html", action="store_true", help="Make html after comparison. Only for parameterization. You must choise --save_fig too.")
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
    elif args.mode == "clusters":
        if args.charges is None or args.sdf_input is None or args.clusters is None or args.atom_type_for_clusterize is None:
            parser.error("For clusterization must be choisen --charges, --sdf_input, --atom_type_for_clusterization and --clusters")
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
    return args, logger
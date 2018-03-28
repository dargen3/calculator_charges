from sys import exit
from termcolor import colored
from os import system
from glob import glob
from modules.parameterization import parameterize


def find_argumets_and_parameterize(path, logger, mode, method_parameterization_args, num_of_parameterized_mol_args, cpu, time):
    par_files = glob("{}*.par".format(path))
    if len(par_files) != 1:
        exit(colored("\n\n\nThere is not 1 parameters file!\n\n\n", "red"))
    sdf_files = glob("{}*.sdf".format(path))
    if len(sdf_files) != 1:
        exit(colored("\n\n\nThere is not 1 sdf file1!\n\n\n", "red"))
    chg_files = glob("{}*.chg".format(path))
    if len(chg_files) != 2:
        exit(colored("\n\n\nThere is not 2 charges files!\n\n\n", "red"))
    sdf_file = sdf_files[0]
    parameters = par_files[0]
    method = parameters.split("_")[-2]
    right_charges = "{}.chg".format(sdf_file[:-4])
    chg_output = "{}_{}.chg".format(sdf_file.split("/")[-1][:-4], method)
    new_parameters = parameters.split("/")[-1]
    num_of_par_mol = None
    if num_of_parameterized_mol_args:
        num_of_par_mol = num_of_parameterized_mol_args
    if method_parameterization_args:
        method_parameterization = method_parameterization_args
    if mode == "parameterization_find_args":
        parameterize(method, parameters, sdf_file, num_of_par_mol, False, right_charges, method_parameterization, new_parameters,
                     chg_output, None, logger, True, True, True, cpu)
    elif mode == "parameterization_send_meta":
        command = "./calculator_charges.py --mode parameterization --method {} --method_parameterization {} --parameters {} --sdf_input {} --right_charges {} " \
                  " --new_parameters {} --chg_output {} -f -v --make_html --save_fig --cpu {} > output.txt 2>&1"\
            .format(method, method_parameterization, parameters.split("/")[-1], sdf_file.split("/")[-1], right_charges.split("/")[-1], new_parameters, chg_output, cpu)
        if num_of_par_mol:
            command += " --num_of_parameterized_mol {}".format(num_of_par_mol)
        system("./modules/send_meta.sh {} {} {} '{}' {} {}".format(parameters, sdf_file, right_charges, command, cpu, time))



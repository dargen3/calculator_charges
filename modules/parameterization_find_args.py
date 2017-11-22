from sys import exit
from termcolor import colored
from os import system
from glob import glob
from modules.parameterization import parameterize


def find_argumets(path, logger):
    par_files = glob("{}*.par".format(path))
    if len(par_files) != 1:
        exit(colored("\n\n\nThere is more parameters file than 1!\n\n\n", "red"))
    sdf_files = glob("{}*.sdf".format(path))
    if len(sdf_files) != 1:
        exit(colored("\n\n\nThere is more sdf file than 1!\n\n\n", "red"))
    sdf_file = sdf_files[0]
    parameters = par_files[0]
    method = parameters.split("_")[1]
    right_charges = "{}.chg".format(sdf_file[:-4])
    chg_output = "{}_{}.chg".format(sdf_file.split("/")[-1][:-4], method)
    new_parameters = parameters.split("/")[-1]
    with open(parameters, "r") as par:
        for line in par:
            try:
                ls = line.split()
                if ls[0]+ls[1]+ls[2]+ls[3] == "Numberofparameterizedmolecules:":
                    num_of_par_mol = int(ls[4])
                    break
            except IndexError:
                pass
    parameterize(method, parameters, sdf_file, num_of_par_mol, False, right_charges, "minimize", new_parameters,
                 chg_output, None, logger, True, True, True, True)



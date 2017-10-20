from modules.set_of_molecule import Set_of_molecule
import os.path
from termcolor import colored
from sys import stdin


def is_the_same(charges, sdf_input, logger):
    if not os.path.isfile(charges):
        exit(colored("There is no charges file with name " + charges + "\n", "red"))
    if not os.path.isfile(sdf_input):
        exit(colored("There is no sdf file with name " + sdf_input + "\n", "red"))
    charges_data = {}
    charges_names = []
    number_of_atoms_chg = 0
    with open(charges, "r+") as charges_file:
        number_of_lines = len(charges_file.readlines())
        number_of_actual_lines = 0
        charges_file.seek(0)
        while number_of_actual_lines < number_of_lines:
            list_with_mol_data = []
            line = charges_file.readline()
            number_of_actual_lines += 1
            name = line.split()[0]
            line = charges_file.readline()
            number_of_actual_lines += 1
            number = int(line.split()[0])
            actual_number = 0
            while actual_number < number:
                line = charges_file.readline()
                number_of_actual_lines += 1
                actual_number += 1
                list_with_mol_data.append(line.split()[1])
                number_of_atoms_chg += 1
            if number_of_actual_lines != number_of_lines:
                charges_file.readline()
                number_of_actual_lines += 1
            charges_data[name] = list_with_mol_data
            charges_names.append(name)
    logger.info(str(charges))
    logger.info("Number of atoms is: " + str(number_of_atoms_chg))
    logger.info("Number of molecules is: " + str(len(charges_data)) + "\n\n\n")
    setm = Set_of_molecule(sdf_input)
    number_of_molecules = len(setm)
    set_of_molecule = setm[:number_of_molecules]
    sdf_data = {}
    sdf_names = []
    number_of_atoms_sdf = 0
    number_of_molecule_sdf = 0
    for molecule in set_of_molecule:
        number_of_atoms_sdf += len(molecule)
        number_of_molecule_sdf += 1
        name_of_molecule = molecule.name
        sdf_names.append(name_of_molecule)
        sdf_data[name_of_molecule] = molecule.atoms_types
    logger.info(sdf_input)
    logger.info("Number of atoms is: " + str(number_of_atoms_sdf))
    logger.info("Number of molecules is: " + str(number_of_molecule_sdf) + "\n\n\n")
    logger.info(colored("Number of together molecules: " + str(len(set(charges_names) & set(sdf_names))) + "\n\n\n",
                        "yellow"))
    return charges_data, sdf_data, charges_names, sdf_names, setm, number_of_lines


def create_new_files(sdf_names, charges_names, sdf_input, charges, setm, number_of_lines):
    together_names = set(sdf_names) & set(charges_names)
    print("Please, write name of new sdf file:")
    new_sdf_file = stdin.readline().rstrip('\n')
    if len(new_sdf_file) == 0:
        exit(colored("File must have same name!", "red"))
    print("Please, write name of new charges file:")
    new_chg_file = stdin.readline().rstrip('\n')
    if len(new_chg_file) == 0:
        exit(colored("File must have same name!", "red"))
    with open(sdf_input, "r") as sdf_input:
        with open(new_sdf_file, "w") as new_sdf_file:
            actual_line = 0
            while actual_line < setm.num_of_lines:
                line = sdf_input.readline()
                actual_line += 1
                if line.split()[0].strip() not in together_names:
                    while line.strip() != "$$$$":
                        line = sdf_input.readline()
                        actual_line += 1
                else:
                    while line.strip() != "$$$$":
                        new_sdf_file.write(line)
                        line = sdf_input.readline()
                        actual_line += 1
                    new_sdf_file.write("$$$$\n")
    with open(charges, "r") as charges:
        actual_line = 0
        charges_dict = {}
        while actual_line < number_of_lines:
            line = charges.readline()
            actual_line += 1
            if line.split()[0] not in together_names:
                while line.strip() != "":
                    line = charges.readline()
                    actual_line += 1
            else:
                list_of_charges = []
                while line.strip() != "":
                    list_of_charges.append(line.strip())
                    line = charges.readline()
                    actual_line += 1
                charges_dict[list_of_charges[0]] = list_of_charges[1:]
    with open(new_chg_file, "w") as new_chg_file:
        for name in sdf_names:
            if name in together_names:
                new_chg_file.write(name + "\n")
                for line in charges_dict[name]:
                    new_chg_file.write(line + "\n")
                new_chg_file.write("\n")

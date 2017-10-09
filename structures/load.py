from termcolor import colored
from sys import exit
from pprint import pprint

def load(file):
    with open(file, 'r') as if_sdf:
        num_of_molecules = 0
        for line in if_sdf:
            if line.strip() == "$$$$":
                num_of_molecules += 1
        if line.strip() != "$$$$":
            exit(colored("ERROR! End of sdf file must be: $$$$", "red"))
    list_with_all_molecules = []
    actual_mol = 0
    with open(file, 'r') as sdf_file:
        while actual_mol < num_of_molecules:
            name = sdf_file.readline().strip()
            for x in range(3):
                line =sdf_file.readline()
            if str(line[-5]) == "2":
                number_of_atoms = int(line[:3])
                number_of_bonds = int(line[3:6])
                list_of_atoms = []
                for x in range(number_of_atoms):
                    line = sdf_file.readline().split()
                    list_of_atoms.append((line[3], x+1, float(line[0]), float(line[1]), float(line[2])))
                list_of_bond = []
                for x in range(number_of_bonds):
                    line = sdf_file.readline().split()
                    if line[2] == "1":
                        type_of_bond = "SINGLE"
                    elif line[2] == "2":
                        type_of_bond = "DOUBLE"
                    elif line[2] == "3":
                        type_of_bond = "TRIPLE"
                    else:
                        exit(colored("ERROR! Wrong type of bond. Multiplicity of bond is " + str(line[2]), "red"))
                    list_of_bond.append((line[0], line[1], type_of_bond))
                line = sdf_file.readline().strip()
                charge = 0
                while line != "$$$$":
                    line = sdf_file.readline().strip()
                    try:
                        l=line.split()
                        if l[0] == "M" and l[1] == "CHG":
                            for x in range(5,3+l[2]*2,2):
                                charge = charge + l[x-1]
                    except IndexError:
                        pass
                dict_with_mol_data = {}
                dict_with_mol_data["atoms"] = list_of_atoms
                dict_with_mol_data["bonds"] = list_of_bond
                dict_with_mol_data["total_charge"] = charge
                list_with_all_molecules.append((name,dict_with_mol_data))
                actual_mol += 1
            elif str(line[-5]) == "3":
                list_of_atoms=[]
                for x in range(3):
                    sdf_file.readline()
                line = sdf_file.readline().split()
                charge=0
                while line[2] != "END" and line[3] != "ATOM":
                    for x in line:
                        if x[:3] == "CHG":
                            charge = charge + float(x[4:])

                    list_of_atoms.append((line[3], int(line[2]), float(line[4]), float(line[5]), float(line[6])))
                    line = sdf_file.readline().split()
                sdf_file.readline()
                line = sdf_file.readline().split()
                list_of_bond = []
                while line[2] != "END" and line[3] != "BOND":
                    if line[3] == "1":
                        type_of_bond = "SINGLE"
                    elif line[3] == "2":
                        type_of_bond = "DOUBLE"
                    elif line[3] == "3":
                        type_of_bond = "TRIPLE"
                    else:
                        exit(colored("ERROR! Wrong type of bond. Multiplicity of bond is " + str(line[3]), "red"))
                    list_of_bond.append((line[5], line[4], type_of_bond))
                    line = sdf_file.readline().split()
                while line != "$$$$":
                    line = sdf_file.readline().strip()
                actual_mol += 1
                dict_with_mol_data = {}
                dict_with_mol_data["atoms"] = list_of_atoms
                dict_with_mol_data["bonds"] = list_of_bond
                dict_with_mol_data["total_charge"] = charge
                list_with_all_molecules.append((name,dict_with_mol_data))


            else:
                print(line)
                exit(colored("ERROR! Sdf file is not correct!", "red"))

    pprint(list_with_all_molecules)



load("proteins/prot_struc.sdf")

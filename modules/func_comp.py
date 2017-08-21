from numpy import sqrt
from scipy import stats
from matplotlib import pyplot as plt
import os.path


def making_dictionary_with_charges(file_with_charges):
    charges_data = {}
    with open(file_with_charges, "r+") as charges:
        number_of_lines = len(charges.readlines())
        number_of_actual_lines = 0
        charges.seek(0)
        while number_of_actual_lines < number_of_lines:
            list_with_mol_data = []
            line = charges.readline()
            number_of_actual_lines += 1
            name = line.split()[0]
            line = charges.readline()
            number_of_actual_lines += 1
            number = int(line.split()[0])
            actual_number = 0
            while actual_number < number:
                line = charges.readline()
                number_of_actual_lines += 1
                actual_number += 1
                list_with_mol_data.append(tuple([line.split()[1], line.split()[2]]))
            if number_of_actual_lines != number_of_lines:
                charges.readline()
                number_of_actual_lines += 1
            charges_data[name] = list_with_mol_data
    return charges_data


def making_final_list(dict_with_charges1, dict_with_charges2):
    dict_with_names1 = dict_with_charges1.keys()
    dict_with_names2 = dict_with_charges2.keys()
    list_with_name = []
    for name in dict_with_names1:
        if name in dict_with_names2:
            list_with_name.append(name)
    list_with_wrong_names = []
    for name in list_with_name:
        for x in range(len(dict_with_charges1[name])):
            if dict_with_charges1[name][x][0] != dict_with_charges2[name][x][0]:
                list_with_wrong_names.append(name)
                break
    for name in list_with_wrong_names:
        list_with_name.remove(name)
        print("Molecule", name, "is not comparised!")
    final_list = []
    list_of_atoms = []
    for name in list_with_name:
        for x in range(len(dict_with_charges1[name])):
            final_list.append((dict_with_charges1[name][x][0], float(dict_with_charges1[name][x][1]), float(dict_with_charges2[name][x][1])))
            list_of_atoms.append(dict_with_charges1[name][x][0])
    final_dict_with_mol = {}
    for name in list_with_name:
        list_with_mol = []
        for x in range(len(dict_with_charges1[name])):
            list_with_mol.append((dict_with_charges1[name][x][0], float(dict_with_charges1[name][x][1]), float(dict_with_charges2[name][x][1])))
        final_dict_with_mol[name] = list_with_mol
    set_of_atoms = set(list_of_atoms)
    return final_list, final_dict_with_mol, tuple(set_of_atoms)

def statistics_for_atom_type(atomic_symbol, list_with_data, atoms):
    list_with_charges1 = []
    list_with_charges2 = []
    try:
        list_with_atomic_data = []
        for atom in list_with_data:
            if atom[0] == atomic_symbol:
                list_with_atomic_data.append(atom)
        list_for_RMSD = []
        for atom in list_with_atomic_data:
            list_for_RMSD.append((atom[1]-atom[2])**2)
        RMSD = sqrt((1.0/len(list_for_RMSD))*sum(list_for_RMSD))
        list_for_deviation = []
        for atom in list_with_atomic_data:
            list_for_deviation.append((atom[1]-atom[2]))
        max_deviation = max(list_for_deviation)
        average_deviation = sum(list_for_deviation)/len(list_for_deviation)
        list_with_charges1 = []
        list_with_charges2 = []
        for atom in list_with_atomic_data:
            list_with_charges1.append(atom[1])
            list_with_charges2.append(atom[2])
        person = stats.pearsonr(list_with_charges1, list_with_charges2)[0]
    except:
        RMSD, max_deviation, average_deviation, person = "-", "-", "-", "-"
    colors = ["#000000", "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#00FFFF", "#FF00FF", "#C0C0C0", "#800000", "#008000", "#800080", "#008080", "#000080"] * 10
    plt.scatter(list_with_charges1, list_with_charges2, marker="o", color=colors[atoms.index(atomic_symbol)], label=atomic_symbol)
    return [atomic_symbol, RMSD, max_deviation, average_deviation, person]


def statistics_for_all_atoms(list_with_data):
    list_for_RMSD = []
    for atom in list_with_data:
        list_for_RMSD.append((atom[1] - atom[2]) ** 2)
    RMSD = sqrt((1.0 / len(list_for_RMSD)) * sum(list_for_RMSD))
    list_for_deviation = []
    for atom in list_with_data:
        list_for_deviation.append((atom[1] - atom[2]))
    max_deviation = max(list_for_deviation)
    average_deviation = sum(list_for_deviation) / len(list_for_deviation)
    list_with_charges1 = []
    list_with_charges2 = []
    for atom in list_with_data:
        list_with_charges1.append(atom[1])
        list_with_charges2.append(atom[2])
    person = stats.pearsonr(list_with_charges1, list_with_charges2)[0]
    return [[RMSD, max_deviation, average_deviation, person]], len(list_with_data)


def average(list):
    return sum(list)/len(list)


def statistics_for_molecules(dictionary):
    list_with_RMSD = []
    list_with_max_deviation = []
    list_with_average_deviation = []
    list_with_pearson = []
    for key in dictionary:
        data = statistics_for_all_atoms(dictionary[key])[0]
        list_with_RMSD.append(data[0][0])
        list_with_max_deviation.append(data[0][1])
        list_with_average_deviation.append(data[0][2])
        list_with_pearson.append(data[0][3])
    return [[average(list_with_RMSD), average(list_with_max_deviation), average(list_with_average_deviation), average(list_with_pearson)]], len(dictionary)


def plotting(charges1, charges2, save_fig):
    plt.xlabel(charges1)
    plt.ylabel(charges2)
    plt.title("Correlation graph")
    plt.xlim([-2, 3])
    plt.ylim([-2, 3])
    x = [-100, 100]
    y = [-100, 100]
    plt.legend()
    plt.plot(x, y)
    if save_fig is not None:
        plt.savefig(save_fig)
    plt.show()


def control_if_arguments_files_exist_for_com(charges1, charges2, save_fig):
    if not os.path.isfile(charges1):
        exit(colored("There is no charges file with name " + charges1 + "\n", "red"))
    if not os.path.isfile(charges2):
        exit(colored("There is no charges file with name " + charges2 + "\n", "red"))
    if save_fig is None:
        if os.path.isfile(str(save_fig) + ".png"):
            print(colored("Warning. There is some file with have the same name like your saved picture from comparison!", "red"))
            print("If you want to replace exist file, please write yes and press enter. Else press enter.")
            decision = stdin.readline().rstrip('\n')
            if decision == "yes":
                os.remove(save_fig + ".png")
                print(colored("Exist file was removed.\n\n\n", "green"))
            else:
                print("\n\n")
                exit(1)

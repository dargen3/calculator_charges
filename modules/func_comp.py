from numpy import sqrt, polyfit, array
from scipy import stats
from matplotlib import pyplot as plt
import os.path
from sys import stdin
from termcolor import colored



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

def making_dictionary_with_charges_para(file_with_charges, atomic_types):
    charges_data = {}
    with open(file_with_charges, "r+") as charges:
        number_of_lines = len(charges.readlines())
        number_of_actual_lines = 0
        charges.seek(0)
        number_of_types = 0
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
                list_with_mol_data.append(tuple([atomic_types[number_of_types], line.split()[2]]))
                number_of_types += 1
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
        print(colored("Molecule", name, "is not comparised!", "red"))
    final_list = []
    list_of_atoms = []
    for name in list_with_name:
        for x in range(len(dict_with_charges1[name])):
            final_list.append((dict_with_charges1[name][x][0], float(dict_with_charges1[name][x][1]),
                               float(dict_with_charges2[name][x][1])))
            list_of_atoms.append(dict_with_charges1[name][x][0])
    final_dict_with_mol = {}
    for name in list_with_name:
        list_with_mol = []
        for x in range(len(dict_with_charges1[name])):
            list_with_mol.append((dict_with_charges1[name][x][0], float(dict_with_charges1[name][x][1]),
                                  float(dict_with_charges2[name][x][1])))
        final_dict_with_mol[name] = list_with_mol
    set_of_atoms = set(list_of_atoms)
    final_dict_with_atoms = {}
    for x in set_of_atoms:
        final_dict_with_atoms[x] = []
    for x in final_list:
        final_dict_with_atoms[x[0]].append(x)
    return final_list, final_dict_with_atoms, final_dict_with_mol, tuple(set_of_atoms)


def statistics_for_atom_type(atomic_symbol, list_with_atomic_data, atoms, fig_all, charges1, charges2, save_fig,
                             axis_range, name_of_all_set):
    try:
        list_for_rmsd = []
        for atom in list_with_atomic_data:
            list_for_rmsd.append((atom[1]-atom[2])**2)
        rmsd = sqrt((1.0/len(list_for_rmsd))*sum(list_for_rmsd))
        list_for_deviation = []
        for atom in list_with_atomic_data:
            list_for_deviation.append(abs(atom[1]-atom[2]))
        max_deviation = max(list_for_deviation)
        average_deviation = sum(list_for_deviation)/len(list_for_deviation)
        list_with_charges1 = []
        list_with_charges2 = []
        for atom in list_with_atomic_data:
            list_with_charges1.append(atom[1])
            list_with_charges2.append(atom[2])
        person_2 = stats.pearsonr(list_with_charges1, list_with_charges2)[0]**2
    except ZeroDivisionError:
        return [atomic_symbol, "-", "-", "-", "-", "-"]
    if save_fig is not False:
        fig_x = plt.figure(atomic_symbol, figsize=(11, 9))
        figx = fig_x.add_subplot(111, rasterized=True)
        figx.set_title(atomic_symbol)
        figx.set_xlabel(charges2, fontsize=15)
        figx.set_ylabel(charges1, fontsize=15)
        figx.set_xlim(axis_range)
        figx.set_ylim(axis_range)
        m, b = polyfit(list_with_charges1, list_with_charges2, 1)
        x = array([x / 20.0 for x in range(-200, 200)])
        figx.plot(x, x, "-", label="y = x", color="black", linewidth=0.5)
        figx.plot(x, m * x + b, ":", color="black", linewidth=0.5, label="y = " + str(round(m, 6))[:4] + "*x +" +
                                                                         str(round(b, 6))[:4])
        figx.plot(list_with_charges1, list_with_charges2, ".", label=atomic_symbol, color="black")
        figx.legend()
        plt.text(axis_range[1], axis_range[0], "Num. of atoms: " + str(len(list_with_charges1)) + "\nrmsd: " +
                 str(rmsd)[:6] + "\nPearson**2: " + str(person_2)[:6], ha='right', va='bottom', fontsize=15)
        name = "{}-{}".format(name_of_all_set, atomic_symbol)
        plt.savefig(name)
    colors = ["#000000", "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#00FFFF", "#FF00FF", "#C0C0C0", "#800000",
              "#008000", "#800080", "#008080", "#000080"] * 10
    fig_all.scatter(list_with_charges1, list_with_charges2, marker=".", color=colors[atoms.index(atomic_symbol)],
                    label=atomic_symbol)
    return [atomic_symbol, rmsd, max_deviation, average_deviation, person_2, len(list_with_atomic_data)]

def statistics_for_all_atoms(list_with_data, fitting=False):
    list_for_rmsd = []
    if fitting:
        min1 = 10
        min2 = 10
        max1 = -10
        max2 = -10
        for atom in list_with_data:
            list_for_rmsd.append((atom[1] - atom[2]) ** 2)
            if atom[1] < min1:
                min1 = atom[1]
            if atom[1] > max1:
                max1 = atom[1]
            if atom[2] < min2:
                min2 = atom[2]
            if atom[2] > max2:
                max2 = atom[2]
        mini = min(min1, min2) - 0.25
        maxi = max(max1, max2) + 0.25
    else:
        for atom in list_with_data:
            list_for_rmsd.append((atom[1] - atom[2]) ** 2)
    rmsd = sqrt((1.0 / len(list_for_rmsd)) * sum(list_for_rmsd))
    list_for_deviation = []
    for atom in list_with_data:
        list_for_deviation.append(abs(atom[1] - atom[2]))
    max_deviation = max(list_for_deviation)
    average_deviation = sum(list_for_deviation) / len(list_for_deviation)
    list_with_charges1 = []
    list_with_charges2 = []
    for atom in list_with_data:
        list_with_charges1.append(atom[1])
        list_with_charges2.append(atom[2])
    if fitting:
        m, b = polyfit(list_with_charges1, list_with_charges2, 1)
        x = array([x/20.0 for x in range(-200, 200)])
        plt.plot(x, x, "-", label="y = x", color="black", linewidth=0.5)
        plt.plot(x, m*x + b, ":", color="black", linewidth=0.5, label="y = " + str(round(m, 6))[:4] + "*x +" +
                                                                      str(round(b, 6))[:4])
    person_2 = stats.pearsonr(list_with_charges1, list_with_charges2)[0]**2
    if fitting:
        return [[rmsd, max_deviation, average_deviation, person_2, len(list_with_data)]], [mini, maxi], \
               len(list_with_charges1)
    return [[rmsd, max_deviation, average_deviation, person_2, len(list_with_data)]]


def average(list_of_number):
    return sum(list_of_number)/len(list_of_number)


def statistics_for_molecules(dictionary, mol_into_log):
    list_with_rmsd = []
    list_with_max_deviation = []
    list_with_average_deviation = []
    list_with_pearson_2 = []
    for key in dictionary:
        data = statistics_for_all_atoms(dictionary[key])
        list_with_rmsd.append(data[0][0])
        list_with_max_deviation.append(data[0][1])
        list_with_average_deviation.append(data[0][2])
        list_with_pearson_2.append(data[0][3])
        if mol_into_log is not None:
            with open(mol_into_log, "a") as mol_log:
                mol_log.write(str(key) + " " + str(data[0][0])[:5] + " " + str(data[0][1])[:5] + " " +
                              str(data[0][2])[:5] + " " + str(data[0][3])[:5] + "\n")
    return [[average(list_with_rmsd), average(list_with_max_deviation), average(list_with_average_deviation),
             average(list_with_pearson_2), len(dictionary)]]


def plotting(charges1, charges2, save_fig_bool, save_fig_name, fig_all, fig, rmsd, person_2, axis_range, number_of_atoms):
    fig_all.set_xlabel(charges2, fontsize=15)
    fig_all.set_ylabel(charges1, fontsize=15)
    fig_all.set_title("Correlation graph", fontsize=15)
    fig_all.set_xlim(axis_range)
    fig_all.set_ylim(axis_range)
    fig_all.legend()
    fig_all.text(axis_range[1], axis_range[0], "Num. of atoms: " + str(number_of_atoms) + "\nrmsd: " + str(rmsd)[:6] +
                 "\nPearson**2: " + str(person_2)[:6], ha='right', va='bottom', fontsize=15)
    if save_fig_bool is not False:
        fig.savefig(save_fig_name)
        plt.show()


def control_if_arguments_files_exist_for_com(charges1, charges2, save_fig, save_fig_true, log, force):
    if not os.path.isfile(charges1):
        exit(colored("There is no charges file with name " + charges1 + "\n", "red"))
    if not os.path.isfile(charges2):
        exit(colored("There is no charges file with name " + charges2 + "\n", "red"))
    try:
        if os.path.isfile(save_fig) and save_fig_true:
            if not force:
                print(colored("Warning. There is some file with have the same name like your saved picture from" +
                              " comparison!", "red"))
                print("If you want to replace exist file, please write yes and press enter. Else press enter.")
                decision = stdin.readline().rstrip('\n')
                if decision == "yes":
                    os.remove(save_fig)
                    print(colored("Exist file was removed.\n\n\n", "green"))
                else:
                    print("\n\n")
                    exit(1)
    except TypeError:
        pass
    try:
        if os.path.isfile(log):
            if not force:
                print(colored("Warning. There is some file with have the same name like your log file from comparison!",
                              "red"))
                print("If you want to replace exist file, please write yes and press enter. Else press enter.")
                decision = stdin.readline().rstrip('\n')
                if decision == "yes":
                    os.remove(log)
                    print(colored("Exist file was removed.\n\n\n", "green"))
                else:
                    print("\n\n")
                    exit(1)
    except TypeError:
        pass

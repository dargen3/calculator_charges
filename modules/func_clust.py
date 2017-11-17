import os.path
from sys import exit
from termcolor import colored
from modules.set_of_molecule import Set_of_molecule
from operator import itemgetter
from matplotlib import pyplot as plt
import seaborn as sns
from numpy import array, unique
from sklearn.cluster import KMeans


def control_if_arguments_files_exist_for_clust(charges, sdf_input):
    if not os.path.isfile(charges):
        exit(colored("There is no parameters file with name " + charges + "\n", "red"))
    if not os.path.isfile(sdf_input):
        exit(colored("There is no sdf file with name " + sdf_input + "\n", "red"))

def check_type_of_atom(atom_type_for_clusterization):
    if len(atom_type_for_clusterization.split("~")) == 1:
        return "atom"
    if len(atom_type_for_clusterization.split("~")) == 2:
        return "atom~high_bond"

def loadig_data_from_sdf_and_chg_file(sdf_input, charges, atom_type, atom_type_for_clusterization):
    locations = {}
    setm = Set_of_molecule(sdf_input)
    number_of_molecules = len(setm)
    set_of_molecule = setm[:number_of_molecules]
    for molecule in set_of_molecule:
        m_symbols = molecule.symbols(atom_type)
        for index, x in enumerate(m_symbols):
            if x == atom_type_for_clusterization:
                if molecule.name not in locations:
                    locations[molecule.name] = []
                bonded_atom_type = []
                for atom in molecule.bonded_atoms[index + 1]:
                    bonded_atom_type.append(m_symbols[atom - 1])
                locations[molecule.name].append((index + 1, bonded_atom_type))
    with open(charges, "r") as charges:
        charges_data = []
        indexes = []
        for line in charges:
            try:
                l0 = line.split()[0]
                if l0 in locations:
                    name = l0
                    indexes = [x[0] for x in locations[l0]]
            except IndexError:
                pass
            l = line.split()
            if len(l) == 3 and int(l[0]) in indexes:
                charges_data.append([name, indexes[0], float(l[2])])
                try:
                    indexes.pop(0)
                except IndexError:
                    pass
    charges_data = sorted(charges_data)
    index = 0
    for x in sorted(locations):
        if len(locations[x]) == 1:
            charges_data[index].append(locations[x][0][1])
            index += 1
        else:
            for y in range(len(locations[x])):
                charges_data[index].append(locations[x][y][1])
                index += 1
    charges_data = sorted(charges_data, key=itemgetter(2))
    charges_values = [x[2] for x in charges_data]
    return charges_data, charges_values

def graph_of_all_charges(atom_type_for_clusterization, charges_values, bw, save_fig):
    if bw is None:
        bw = 0.1
    else:
        bw = float(bw)
    fig = plt.figure(figsize=(11, 9))
    fig_s = fig.add_subplot(111)
    title = "Kernel density plot for {}".format(atom_type_for_clusterization)
    fig_s.set_title(title, fontsize=15)
    fig_s.set_xlabel("Charge", fontsize=12)
    fig_s.set_ylabel("Density (bw = {})".format(bw), fontsize=12)
    sns.set_style('whitegrid')
    sns.kdeplot(array(charges_values), bw=bw)
    if save_fig:
        fig.savefig(title.replace(" ", "_"))
    plt.show()
    exit()

def clusterize(charges_values, clusters, charges_data):
    charges_values_2D = [[x, x] for x in charges_values]
    kmeans = KMeans(n_clusters=int(clusters), random_state=0).fit(charges_values_2D)
    distribution = kmeans.labels_
    items, count = unique(distribution, return_counts=True)
    counts = dict(zip(items, count))
    for index, x in enumerate(distribution):
        charges_data[index].append(x)
    average_charges = kmeans.cluster_centers_
    return counts, average_charges, charges_data

def print_clusters(clusters, counts, average_charges, charges_data):
    for x in range(int(clusters)):
        print("Cluster: {}".format(x + 1))
        num_of_atoms = counts[x]
        print("Number of atoms: {}".format(num_of_atoms))
        print("Average charge {:.2f}:".format(average_charges[x][0]))
        atoms_count = {}
        for charge in charges_data:
            if charge[4] == x:
                for bonded_atom in charge[3]:
                    if bonded_atom not in atoms_count:
                        atoms_count[bonded_atom] = 1
                    else:
                        atoms_count[bonded_atom] += 1
        stat = []
        for x in atoms_count:
            stat.append((x, float(atoms_count[x] / num_of_atoms)))
        stat = sorted(stat, key=itemgetter(1))[::-1]
        average_num_of_bonded_atoms = 0
        for x in stat:
            average_num_of_bonded_atoms += x[1]
        print("Average number of bonded atoms: {:.2f}".format(average_num_of_bonded_atoms))
        for index, x in enumerate(stat):
            if index == 0:
                print("Average atom composition: ", end="")
                print(" {}: {:.2f}".format(x[0], float(x[1])))
            else:
                print("{:>28}: {:.2f}".format(x[0], float(x[1])))
        print("\n\n\n\n")
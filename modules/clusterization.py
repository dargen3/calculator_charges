import os.path
from sys import exit
from termcolor import colored
from .set_of_molecule import Set_of_molecule
from operator import itemgetter
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde
from numpy import array, unique, linspace
from sklearn.cluster import KMeans
from .statistics import statistics


def control_if_arguments_files_exist_for_clust(charges, sdf_input):
    if not os.path.isfile(charges):
        exit(colored("There is no charges file with name " + charges + "\n", "red"))
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
    print(len(locations))
    switch = True
    with open(charges, "r") as charges:
        charges_data = []
        indexes = []
        for line in charges:
            try:
                l0 = line.split()[0]
                if l0 in locations and switch:
                    name = l0
                    indexes = [x[0] for x in locations[l0]]
            except IndexError:
                switch = True
                continue
            ls = line.split()
            if len(ls) == 3 and int(ls[0]) in indexes:
                charges_data.append([name, indexes[0], float(ls[2])])
                try:
                    indexes.pop(0)
                except IndexError:
                    pass
            switch = False
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
    kde = gaussian_kde(charges_values, bw_method=bw)
    x_grid = linspace(min(charges_values) - (max(charges_values) - min(charges_values)) * 0.3,
                      max(charges_values) + (max(charges_values) - min(charges_values)) * 0.3, 1000)
    kde.evaluate(x_grid)
    fig_s.plot(x_grid, kde(x_grid))
    if save_fig:
        fig.savefig(title.replace(" ", "_"))
    plt.show()
    exit()


def clust(charges_values, clusters, charges_data):
    charges_values_2d = [[x, x] for x in charges_values]
    kmeans = KMeans(n_clusters=int(clusters), random_state=0).fit(charges_values_2d)
    distribution = kmeans.labels_
    items, count = unique(distribution, return_counts=True)
    counts = dict(zip(items, count))
    for index, x in enumerate(distribution):
        charges_data[index].append(x)
    average_charges = kmeans.cluster_centers_
    return counts, average_charges, charges_data


def print_clusters(clusters, counts, average_charges, charges_data):
    for cluster in range(int(clusters)):
        print("Cluster: {}".format(cluster + 1))
        num_of_atoms = counts[cluster]
        print("Number of atoms: {}".format(num_of_atoms))
        print("Average charge {:.2f}:".format(average_charges[cluster][0]))
        atoms_count = {}
        for charge in charges_data:
            if charge[4] == cluster:
                for bonded_atom in charge[3]:
                    if bonded_atom not in atoms_count:
                        atoms_count[bonded_atom] = 1
                    else:
                        atoms_count[bonded_atom] += 1
        stat = []
        for number in atoms_count:
            stat.append((number, float(atoms_count[number] / num_of_atoms)))
        stat = sorted(stat, key=itemgetter(1))[::-1]
        average_num_of_bonded_atoms = 0
        for x in stat:
            average_num_of_bonded_atoms += x[1]
        print("Average number of bonded atoms: {:.2f}".format(average_num_of_bonded_atoms))
        for index, comp in enumerate(stat):
            if index == 0:
                print("Average atom composition: ", end="")
                print(" {}: {:.2f}".format(comp[0], float(comp[1])))
            else:
                print("{:>30}: {:.2f}".format(comp[0], float(comp[1])))
        print("\n\n\n\n")


def clusterize(charges, sdf_input, logger, atom_type_for_clusterization, clusters, fine_of_graph, save_fig):
    control_if_arguments_files_exist_for_clust(charges, sdf_input)
    logger.info("Checking, if sdf file and chg file are for the same molecules...\n\n\n")
    statistics(charges, sdf_input, logger)
    logger.info("\n\n\nLoading data... from {} and {}.".format(sdf_input, charges))
    atom_type = check_type_of_atom(atom_type_for_clusterization)
    charges_data, charges_values = loadig_data_from_sdf_and_chg_file(sdf_input, charges, atom_type,
                                                                     atom_type_for_clusterization)

    logger.info(colored("Loading of data was sucessfull. \n\n\n", "green"))
    if clusters == "0":
        graph_of_all_charges(atom_type_for_clusterization, charges_values, fine_of_graph, save_fig)
    logger.info("Clusterization of {} clusters...".format(clusters))
    counts, average_charges, charges_data = clust(charges_values, clusters, charges_data)
    logger.info(colored("Clusterization was sucessfull. \n\n\n", "green"))
    print_clusters(clusters, counts, average_charges, charges_data)

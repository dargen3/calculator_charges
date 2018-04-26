from glob import glob
from os import system
from termcolor import colored
from modules.set_of_molecule import Set_of_molecule


def make_html(name, sdf, method, all_atoms, all_molecules, fig_all, atomic_types, atomic_types_data, time, method_para,
              num_of_par_mol, validation, time_of_parameterization):
    name = name + ".html"
    with open(name, "w") as html_file:
        print("\n\ntest\n\n")
        print(method, sdf, time)
        print("\n\ntest\n\n")
        html_file.write("<h1>Method: {0}</h1>\n".format(method))
        html_file.write("<h1>Set of molecules: {0}</h1>\n".format(sdf))
        html_file.write("<h1>Date of parameterization: {0}</h1>\n".format(time))
        html_file.write("<h1>Time of parameterization: {0}</h1>\n".format(time_of_parameterization))
        html_file.write("<h1>Method of parameterization: {0}</h1>\n".format(method_para))
        html_file.write("<h1>Number of parameterized molecules: {0}</h1>\n".format(num_of_par_mol))
        if validation:
            html_file.write("<h1>Mode: Validation 70:30</h1>\n")
        else:
            html_file.write("<h1>Mode: Full set parameterization</h1>\n")
        html_file.write("<h2>Atoms:</h2>\n")
        html_file.write("<table border=1>\n")
        html_file.write("<tbody>\n")
        html_file.write("<th>Number of atoms</th>\n")
        html_file.write("<th>RMSD</th>\n")
        html_file.write("<th>Pearson**2</th>\n")
        html_file.write("<th>Maximum deviation</th>\n")
        html_file.write("<th>Average deviation</th>\n")
        html_file.write("</tr>\n")
        html_file.write("<tr>\n")
        html_file.write("<td>{0}</td>\n".format(str(all_atoms[4])))
        html_file.write("<td>{0}</td>\n".format(str(all_atoms[0])[:6]))
        html_file.write("<td>{0}</td>\n".format(str(all_atoms[3])[:6]))
        html_file.write("<td>{0}</td>\n".format(str(all_atoms[1])[:6]))
        html_file.write("<td>{0}</td>\n".format(str(all_atoms[2])[:6]))
        html_file.write("</tr>\n")
        html_file.write("<tbody>\n")
        html_file.write("</table>\n")
        html_file.write("<h2>Molecules:</h2>\n")
        html_file.write("<table border=1>\n")
        html_file.write("<tbody>\n")
        html_file.write("<th>Number of molec.</th>\n")
        html_file.write("<th>RMSD</th>\n")
        html_file.write("<th>Pearson**2</th>\n")
        html_file.write("<th>Maximum deviation</th>\n")
        html_file.write("<th>Average deviation</th>\n")
        html_file.write("</tr>\n")
        html_file.write("<tr>\n")
        html_file.write("<td>{0}</td>\n".format(str(all_molecules[4])))
        html_file.write("<td>{0}</td>\n".format(str(all_molecules[0])[:6]))
        html_file.write("<td>{0}</td>\n".format(str(all_molecules[3])[:6]))
        html_file.write("<td>{0}</td>\n".format(str(all_molecules[1])[:6]))
        html_file.write("<td>{0}</td>\n".format(str(all_molecules[2])[:6]))
        html_file.write("</tr>\n")
        html_file.write("<tbody>\n")
        html_file.write("</table>\n")
        # picture, all atoms
        html_file.write("<img src=\"{0}\" width=\"1600\">\n".format(fig_all))
        # table with statistics for atomic types
        html_file.write("<h2>Atomic types:</h2>\n")
        html_file.write("<table border=1>\n")
        html_file.write("<tbody>\n")
        html_file.write("<th>Atomic type</th>\n")
        html_file.write("<th>Number of atoms</th>\n")
        html_file.write("<th>RMSD</th>\n")
        html_file.write("<th>Pearson**2</th>\n")
        html_file.write("<th>Maximum deviation</th>\n")
        html_file.write("<th>Average deviation</th>\n")
        html_file.write("</tr>\n")
        for x in atomic_types_data:
            if x[1] < 0.05:
                html_file.write("<tr style=\"background-color: green;\">\n")
            elif x[1] < 0.1:
                html_file.write("<tr style=\"background-color: #4ca64c;\">\n")
            elif x[1] < 0.15:
                html_file.write("<tr style=\"background-color: #99cc99;\">\n")
            elif x[1] < 0.2:
                html_file.write("<tr style=\"background-color: yellow;\">\n")
            elif x[1] < 0.3:
                html_file.write("<tr style=\"background-color: orange;\">\n")
            elif x[1] < 0.4:
                html_file.write("<tr style=\"background-color: red; color: white;\">\n")
            elif x[1] >= 0.4:
                html_file.write("<tr style=\"background-color: darkred; color: white;\">\n")
            html_file.write("<td>{0}</td>\n".format(x[0]))
            html_file.write("<td>{0}</td>\n".format(str(x[5])))
            html_file.write("<td>{0}</td>\n".format(str(x[1])[:6]))
            html_file.write("<td>{0}</td>\n".format(str(x[4])[:6]))
            html_file.write("<td>{0}</td>\n".format(str(x[2])[:6]))
            html_file.write("<td>{0}</td>\n".format(str(x[3])[:6]))
            html_file.write("</tr>\n")
        html_file.write("<tbody>\n")
        html_file.write("</table>\n")
        html_file.write("<table border=\"1\" style=\"margin-top: 0.5em\">\n")
        html_file.write("<tr>\n")
        html_file.write("<td><strong>Legend:</strong> RMSD </td>\n")
        html_file.write("<td style=\"background-color: green; color: white;\">&lt; 0.05</td>\n")
        html_file.write("<td style=\"background-color: #4ca64c; \">&lt; 0.1</td>\n")
        html_file.write("<td style=\"background-color: #99cc99;\">&lt; 0.15</td>\n")
        html_file.write("<td style=\"background-color: yellow;\">&lt; 0.2</td>\n")
        html_file.write("<td style=\"background-color: orange;\">&lt; 0.3</td>\n")
        html_file.write("<td style=\"background-color: red; color: white;\">&lt; 0.4</td>\n")
        html_file.write("<td style=\"background-color: darkred; color: white;\">&gt;= 0.4</td>\n")
        html_file.write("</tr>\n")
        html_file.write("</table>\n")
        # pictures for types of atoms
        html_file.write("<div style=\"width: 1600px;\">")
        for x in atomic_types:
            path = "{}-{}.png".format(fig_all[:-4], x)
            html_file.write("<img src=\"{}\" style=\"float: left; width: 800px;\">".format(path))


def make_complete_html(verbose):
    if verbose:
        print("Creating html...")
    with open("data/index.html", "w") as html_file:
        html_file.write("<h1>Calculator charges</h1>\n")
        html_file.write("<h2>Source code: https://github.com/dargen3/calculator_charges</h2>\n")
        html_file.write("<a href = \"data\">All data</a>\n<br />\n<br />\n<br />\n<br />\n")
        html_files = glob("data/*/*/*.html")
        methods_data_html = {}
        sets_of_molecules = []
        for file in html_files:
            splitted_path = file.split("/")
            method = splitted_path[1]
            sets_of_molecules.append(splitted_path[-2])
            if method not in methods_data_html:
                methods_data_html[method] = []
            methods_data_html[method].append(splitted_path)
        sdf_files = glob("data/*/*/*.sdf")
        sdf_files_check = []
        for sdf_file_path in sdf_files:
            sdf_file = sdf_file_path.split("/")[-1]
            if sdf_file not in sdf_files_check:
                Set_of_molecule(sdf_file_path).statistics_data(write_to_file="data/sets_of_molecules_info/{}_info.txt".format(sdf_file[:-4]))
                sdf_files_check.append(sdf_file)
        sets_int = []
        sets_str = []
        for setm in set(sets_of_molecules):
            try:
                sets_int.append(int(setm))
            except ValueError:
                sets_str.append(setm)
        sets_of_molecules = [str(sets) for sets in sorted([seti for seti in sets_int])] + sorted(sets_str)
        html_file.write("<table border=1>\n")
        html_file.write("<tbody>\n")
        html_file.write("<th>Method</th>\n")
        for setm in sets_of_molecules:
            html_file.write("<th>{}</th>\n".format(setm))
        html_file.write("</tr>\n")
        names = {"1956": "DTP_small", "4475": "DTP_large", "proteins": "", "8144": "CCD_gen_CHNO", "17769": "CCD_gen_all"}
        html_file.write("<th>Name_of_set</th>\n")
        for setm in sets_of_molecules:
            try:
                html_file.write("<th>{}</th>\n".format(names[setm]))
            except KeyError:
                html_file.write("<th></th>\n")
        html_file.write("</tr>\n")
        html_file.write("<tr>\n")
        html_file.write("<td>Set of molecules info</td>")
        for set_info in sets_of_molecules:
            html_file.write("<td><a href = \"data/sets_of_molecules_info/{}_info.txt\">{}_info</a>\n<br />\n</td>\n".format(set_info , set_info))
        html_file.write("</tr>\n")
        for method in sorted(methods_data_html):
            html_file.write("<tr>\n")
            html_file.write("<td>{}</td>\n".format(method))
            for setm in sets_of_molecules:
                is_there = False
                for setm_of in methods_data_html[method]:
                    if setm_of[2] == setm:
                        is_there = True
                        break
                if is_there == True:
                    path = "data/{}/{}/{}_{}.html".format(method, setm, setm, method)
                    name = "{} {}".format(setm, method)
                    html_file.write("<td><a href = \"{}\">{}</a>\n<br />\n</td>\n".format(path, name))
                else:
                    html_file.write("<td>no results</td>\n".format(setm))
            html_file.write("</tr>\n")
        html_file.write("<tbody>\n")
        html_file.write("</table>\n")
        html_file.write("<br /><br /><br /><br /><h3>Contact: dargen3@centrum.cz</h3>\n")
        print("Copying of data...\n\n\n")
    system("scp -r  data dargen3@lcc.ncbr.muni.cz:/home/dargen3/www/")
    if verbose:
        print(colored("\n\n\nData was copied sucessfully.\n\n\n", "green"))
        print("Setting permissions...")
    system("ssh dargen3@lcc.ncbr.muni.cz \" mv www/data/index.html www/index.html ; chmod -R 705 * \"")
    if verbose:
        print(colored("Setting of permissions was sucessfull.\n\n\n", "green"))



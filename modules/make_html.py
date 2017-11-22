from glob import glob
from os import system
from termcolor import colored


def make_html(name, sdf, method, all_atoms, all_molecules, fig_all, atomic_types, atomic_types_data, time, method_para,
              num_of_par_mol, validation):
    name = name + ".html"
    with open(name, "w") as html_file:
        html_file.write("<h1>Method: {0}</h1>\n".format(method))
        html_file.write("<h1>Set of molecules: {0}</h1>\n".format(sdf))
        html_file.write("<h1>Date of parameterization: {0}</h1>\n".format(time))
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
                html_file.write("<tr style=\"background-color: green; color: white;\">\n")
            elif x[1] < 0.1:
                html_file.write("<tr style=\"background-color: #4ca64c; color: white;\">\n")
            elif x[1] < 0.15:
                html_file.write("<tr style=\"background-color: #99cc99; color: white;\">\n")
            elif x[1] < 0.2:
                html_file.write("<tr style=\"background-color: yellow; color: white;\">\n")
            elif x[1] < 0.3:
                html_file.write("<tr style=\"background-color: orange; color: white;\">\n")
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
        html_file.write("<td style=\"background-color: red;\">&lt; 0.4</td>\n")
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
        html_file.write("<h1>Charges method data</h1>\n")
        html_file.write("<a href = \"data\">All data</a>\n<br />\n")
        files = glob("data/*/*/*.html")
        methods = {}
        for file in files:
            splitted_path = file.split("/")
            method = splitted_path[1]
            if method not in methods:
                methods[method] = []
            methods[method].append(splitted_path)
        for method in methods:
            html_file.write("<h2>{}</h2>\n".format(method))
            for result in methods[method]:
                name = result[3][:-5]
                file = ""
                for x in result:
                    file += "/" + x
                print(name)
                html_file.write("<a href = \"{}\">{}</a>\n<br />\n".format(file[1:], name))
    if verbose:
        print(colored("Html was created sucessfully!\n\n\n", "green"))
        print("Copying of data...\n\n\n")
    system("scp -r  data dargen3@lcc.ncbr.muni.cz:/home/dargen3/www/")
    if verbose:
        print(colored("\n\n\nData was copied sucessfully.\n\n\n", "green"))
        print("Setting permissions...")
    system("scp -r  data/index.html dargen3@lcc.ncbr.muni.cz:/home/dargen3/www/")
    system("ssh dargen3@lcc.ncbr.muni.cz \" mv www/data/index.html www/index.html ; chmod -R 705 * \"")
    if verbose:
        print(colored("Setting of permissions was sucessfull.\n\n\n", "green"))



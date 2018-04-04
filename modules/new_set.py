from os import system
from os import path as ospath

def new_set(path, right_charges, sdf_input, parameters):
    if ospath.isfile(right_charges) and ospath.isfile(sdf_input):
        data = path.split("/")
        method = data[1]
        name = data[2]
        parameters_name = "par_{}_{}.par".format(method, name)
        charges = right_charges.replace(".", "_{}.".format(method))
        system("mkdir {}".format(path))
        system("cp {} {}".format(right_charges, charges))
        system("mv {} {} {} {}".format(charges, sdf_input, right_charges, path))
        system("cp {} {}".format(parameters, "{}/{}".format(path, parameters_name)))
        print("New set was created successful!")
        print("Path is {}.".format(path))
    else:
        print("Error! There is no right_charges file or sdf_input file!!")

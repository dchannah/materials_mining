#!/usr/bin/python

import sys
import ast
import math
import yaml
import urllib
import argparse
import numpy as np
import matplotlib.pyplot as plt
from diffusion_updater_v3 import UpdaterNEB
from pymongo import MongoClient
from pymatgen import MPRester, Structure, Site
from scipy import interpolate

# Arbitrary pair index for the diffusion updater
pair_index = "001"

# Label for combined properties
combo_label = "force from closest cat + anion"

# Some global dictionaries to simplify labeling.
label_dict = {"path_energy": "Rel. Energy (meV)", "effective_coordination": "Effective Coord. No.",
              "average_bond_length": "Average Bond Length (A)", "distortion_index": "Distortion Index",
              "coordination": "Coordination Number", "distance": "Distance along path (\%)",
              "activation_energy": "Activation Energy (meV)", "path_length": "Path Length (A)",
              "volume_per_anion": "Volume per O$^{2-}$ (\\r{A}$^3$)", "min_cation_distance": "Cation Distance (\\r{A})",
              "av_delta_cation": "Total CN change", "site_e_diff": "E$_i$ - E$_s$ (meV)",
              "site_cn_diff": "CN$_i$ - CN$_s$", "min_cat_dist": "Nearest cation distance (\\r{A})",
              "combo_prop": combo_label, "min_anion_dist": "Nearest anion distance (\\r{A})",
              "min_anion_distance": "Anion Distance (\\r{A})", "site_index": "Site Index",
              "displacements_by_site": "Total Displacement (\\r{A})"}

color_dict = {"Mg": "red", "Ca": "green", "Zn": "blue", "Li": "orange", "Na": "darkmagenta"}
tm_color_dict = {"Ti": "black", "V": "red", "Cr": "green", "Mn": "violet", "Fe": "white", "Co": "blue",
                 "Ni": "orange", "Mo": "brown"}

# Some global settings to make plots look pretty
plt.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]
plotting_dpi = 100


# Routines
def initialize_database(db_config_file):
    """ Connect to the mongoDB database and pull a collection to work with. """
    with open(db_config_file, 'r') as f:
        db_config = yaml.load(f)
    client = MongoClient("mongodb://" + db_config['user'] + ":" + urllib.quote_plus(db_config['passwd']) + "@" +
                         db_config['host'] + ":" + str(db_config['port']) + "/" + db_config['db'])
    db = client[db_config['db']]
    collection = db[db_config['collection']]
    return collection


def plot_individual(xdict, ydict, xprop, yprop, documents, spline):
    """ Plot two properties from a set of NEB paths against each.
    The most common example is energy vs. distance along the path."""
    figure_array = {}
    for item in documents:
        xlabel = "\\textbf{" + label_dict[xprop] + "}"
        ylabel = "\\textbf{" + label_dict[yprop] + "}"
        x = xdict[item["path_id"]]
        y = ydict[item["path_id"]]
        # fig_title = item["path_id"] + "(" + item["pretty_formula"] + ")"  # Individual traces
        # fig_title = yprop + item["cation_type"]  # Plot by cation
        fig_title = yprop  # All together
        figure_array[item["path_id"]] = plt.figure(fig_title, figsize=(6,6), dpi=plotting_dpi)
        ax = figure_array[item["path_id"]].add_subplot(111) 
        ax.scatter(x,y, s=70, zorder=2, color=color_dict[item["cation_type"]], linewidths=2.5, edgecolors='black')
        if spline:
            tck = interpolate.splrep(x, y, s=0)
            xnew = np.arange(0, 100, 0.1)
            splfit = interpolate.splev(xnew, tck, der=0)
            x = xnew
            y = splfit
        if item["path_id"][-3:] == "002":
            ax.plot(x,y, linewidth=2.5, zorder=1, color=color_dict[item["cation_type"]], linestyle='dashed')
        elif item["path_id"][-3:] == "003":
            ax.plot(x,y, linewidth=2.5, zorder=1, color=color_dict[item["cation_type"]], linestyle='dotted')
        else:
            ax.plot(x,y, linewidth=2.5, zorder=1, color=color_dict[item["cation_type"]])
        ax.set_xlabel(xlabel, fontsize=24)
        # ax.set_ylim([0,1200])
        # ax.set_xlim([0,100])
        ax.set_ylabel(ylabel, fontsize=24)
        ax.tick_params(axis='x', labelsize=22)
        ax.tick_params(axis='y', labelsize=22)
        border_width = 2
        [i.set_linewidth(border_width) for i in ax.spines.itervalues()]
        plt.tight_layout()
        plt.legend(loc='best', prop={'size': 14})
        plt.rc('text', usetex=True)
        plt.rc('font', family='sans-serif')
        plt.tight_layout()
    plt.show()


def plot_individual_tm(xdict, ydict, xprop, yprop, documents, spline):
    """ Plot two properties from a set of NEB paths against each.
    The most common example is energy vs. distance along the path."""
    figure_array = {}
    for item in documents:
        xlabel = "\\textbf{" + label_dict[xprop] + "}"
        ylabel = "\\textbf{" + label_dict[yprop] + "}"
        print str(item["path_id"])
        x = xdict[item["path_id"]]
        y = ydict[item["path_id"]]
        # fig_title = item["path_id"] + "(" + item["pretty_formula"] + ")"  # Individual traces
        fig_title = yprop + item["cation_type"]  # Plot by cation
        figure_array[item["path_id"]] = plt.figure(fig_title, figsize=(6,6), dpi=plotting_dpi)
        ax = figure_array[item["path_id"]].add_subplot(111)
        ax.scatter(x,y, s=70, zorder=2, color=tm_color_dict[item["tm_type"][0]], linewidths=2.5, edgecolors='black',
                   label=item["tm_type"][0])
        if spline:
            tck = interpolate.splrep(x, y, s=0)
            xnew = np.arange(0, 100, 0.1)
            splfit = interpolate.splev(xnew, tck, der=0)
            x = xnew
            y = splfit
        if item["path_id"][-3:] == "002":
            ax.plot(x, y, linewidth=2.5, zorder=1, color=tm_color_dict[item["tm_type"][0]], linestyle='dashed')
        else:
            ax.plot(x, y, linewidth=2.5, zorder=1, color=tm_color_dict[item["tm_type"][0]])
        ax.set_xlabel(xlabel, fontsize=24)
        # ax.set_ylim([0,1200])
        # ax.set_xlim([7,22])
        ax.set_ylabel(ylabel, fontsize=24)
        ax.tick_params(axis='x', labelsize=22)
        ax.tick_params(axis='y', labelsize=22)
        border_width = 2
        [i.set_linewidth(border_width) for i in ax.spines.itervalues()]
        plt.tight_layout()
        plt.legend(loc='best', prop={'size': 14})
        plt.rc('text', usetex=True)
        plt.rc('font', family='sans-serif')
        plt.tight_layout()
    plt.show()


def write_outfile_individual(xdict, ydict, xprop, yprop, documents):
    for item in documents:
        filename = item["path_id"] + "_" + xprop + "_" + yprop + ".dat"
        outfile = open(filename, 'w')
        x_data = xdict[item["path_id"]]
        y_data = ydict[item["path_id"]]
        for i in range(len(x_data)):
            outfile.write(str(x_data[i]) + " " + str(y_data[i]) + "\n")
    return


def plot_collective(xdict, ydict, xprop, yprop, documents):
    """ Plot a single valued property (like barrier height) against another for a range of compounds."""
    x_ion = {"Mg": [], "Ca": [], "Zn": [], "Li": [], "Na": []}
    y_ion = {"Mg": [], "Ca": [], "Zn": [], "Li": [], "Na": []}
    for item in documents:
        if item["path_id"][-3:] == "001":
            x_ion[item["cation_type"]].append(xdict[item["path_id"]])
            y_ion[item["cation_type"]].append(ydict[item["path_id"]])
    fig = plt.figure(figsize=(6,6), dpi=plotting_dpi)
    ax = fig.add_subplot(111)
    for ion in ["Mg", "Ca", "Zn", "Li", "Na"]:
        ax.scatter(x_ion[ion], y_ion[ion], s=70, zorder=2, color=color_dict[ion], linewidths=2.5, edgecolors='black',
                   label=ion)
    xlabel = "\\textbf{" + label_dict[xprop] + "}"
    ylabel = "\\textbf{" + label_dict[yprop] + "}"
    
    # # Plot lines for fitting, if useful
    # x2 = np.arange(-700, 3300, 1)
    # ax.plot(x2, x2)
    
    # # For setting axis boundaries
    # ax.set_xlim([-700, 3500])
    # ax.set_ylim([0,100])
    
    # Plot display settings
    ax.set_xlabel(xlabel, fontsize=24)
    ax.set_ylabel(ylabel, fontsize=24)
    ax.tick_params(axis='x', labelsize=22)
    ax.tick_params(axis='y', labelsize=22)
    border_width = 2
    [i.set_linewidth(border_width) for i in ax.spines.itervalues()]
    plt.tight_layout()
    plt.legend(loc='best', prop={'size':10})
    # plt.legend(loc='best')
    plt.rc('text', usetex=True)
    plt.rc('font', family='sans-serif')
    plt.show()


def plot_collective_tm(xdict, ydict, xprop, yprop, documents):
    """ Plot a single valued property (like barrier height) against another for a range of compounds."""
    x_ion = {"Ti": [], "V": [], "Fe": [], "Cr": [], "Mn": [], "Co": [], "Ni": [], "Mo": []}
    y_ion = {"Ti": [], "V": [], "Fe": [], "Cr": [], "Mn": [], "Co": [], "Ni": [], "Mo": []}
    for item in documents:
        x_ion[item["tm_type"][0]].append(xdict[item["path_id"]])
        y_ion[item["tm_type"][0]].append(ydict[item["path_id"]])
    fig = plt.figure(figsize=(6,6), dpi=plotting_dpi)
    ax = fig.add_subplot(111)
    for ion in ["Ti", "V", "Cr", "Mn", "Co", "Ni"]:
        ax.scatter(x_ion[ion], y_ion[ion], s=70, zorder=2, color=tm_color_dict[ion], linewidths=2.5, edgecolors='black',
                   label=ion)
    xlabel = "\\textbf{" + label_dict[xprop] + "}"
    ylabel = "\\textbf{" + label_dict[yprop] + "}"

    # # Plot lines for fitting, if useful
    # x2 = np.arange(-700, 3300, 1)
    # ax.plot(x2, x2)

    # # For setting axis boundaries
    # ax.set_xlim([-700, 3500])
    # ax.set_ylim([0,100])

    # Plot display settings
    ax.set_xlabel(xlabel, fontsize=24)
    ax.set_ylabel(ylabel, fontsize=24)
    ax.tick_params(axis='x', labelsize=22)
    ax.tick_params(axis='y', labelsize=22)
    border_width = 2
    [i.set_linewidth(border_width) for i in ax.spines.itervalues()]
    plt.tight_layout()
    plt.legend(loc='best', prop={'size':10})
    plt.rc('text', usetex=True)
    plt.rc('font', family='sans-serif')
    plt.show()


def get_quantity(document, property):
    working_ion = document["cation_type"]
    structures_as_dict = document["neb_images"] 
    structures = []
    for i in range(len(structures_as_dict)):
        structures.append(Structure.from_dict(structures_as_dict[i]))
    update = UpdaterNEB(document["mp_id"], document["cation_type"], pair_index)
    diffuser_index = update.get_diffuser_index(structures, working_ion)
    if property in ["path_energy", "path_length", "coordination", "effective_coordination", "average_bond_length",
                    "distortion_index", "activation_energy"]:
        property_value = document["NEB_analysis"][property]
    elif property == "distance":
        property_value = path_percentage(structures, document, diffuser_index)
    elif property == "volume_per_anion":
        structure = structures[0] 
        property_value = vol_per_anion(structure)
    elif property == "min_cation_distance":
        # structure = structures[0]  # For the starting site
        # Below is the code to get the cation distance in the intermediate site.
        neb_energies = property_value = document["NEB_analysis"]["path_energy"]
        activated_state = get_max_index(neb_energies)
        structure = structures[activated_state]
        property_value = get_min_cation_distance(structure, diffuser_index)
    elif property == "min_anion_distance":
        structure_initial = structures[0]  # For the starting site
        # Below is the code to get the cation distance in the intermediate site.
        neb_energies = document["NEB_analysis"]["path_energy"]
        activated_state = get_max_index(neb_energies)
        structure_activated = structures[activated_state]
        property_value = get_min_anion_distance(structure_activated, diffuser_index)
    elif property == "av_delta_cation":
        cn_list = document["NEB_analysis"]["effective_coordination"]
        property_value = calculate_average_change(cn_list)
    elif property == "site_e_diff":
        path_energies = document["NEB_analysis"]["path_energy"]
        property_value = site_energy_difference(path_energies)
        property_value = math.sqrt(property_value**2)
    elif property == "site_cn_diff":
        path_cn = document["NEB_analysis"]["effective_coordination"]
        property_value = site_energy_difference(path_cn)
        property_value = math.sqrt(property_value**2)
    elif property == "min_cat_dist":
        property_value = find_closest_cation(structures, diffuser_index)
    elif property == "min_anion_dist":
        property_value = find_closest_anion(structures, diffuser_index)
    elif property == "site_index":
        sample_struct = structures[0]
        property_value = np.arange(0, len(sample_struct.sites))
    elif property == "displacements_by_site":
        property_value = get_all_displacements(structures)
    elif property == "combo_prop":
        prop1 = get_quantity(document, "min_anion_dist")
        prop2 = get_quantity(document, "min_cat_dist")
        property_value = return_function(prop1, prop2)
    else:
        property_value = document[property]
    return property_value


def path_percentage(structures, document, diffuser_index): 
    path_distance = 0
    path_sum = get_quantity(document, "path_length")
    percentage_list = [0]
    for i in range(0, len(structures)-1):
        path_distance = path_distance + structures[i][diffuser_index].distance(structures[i+1][diffuser_index])
        percentage_list.append(path_distance/path_sum*100)
    return percentage_list


def get_site_displacement(structure_array, index):
    """
    :param structure_array: NEB images array or array along a magnesiation path (for example)
    :param index: Index of the site whose displacement we will sum up along the path
    :return: A sum total of the distance traveled by that site along the path
    """
    distance = 0
    for i in range(len(structure_array)-1):
        current = structure_array[i][index]
        next = structure_array[i+1][index]
        distance += current.distance(next)
    return distance


def get_all_displacements(structure_array):
    """
    :param structure_array: Array of structures representing some sequence (NEB, magnesiation, etc.)
    :return: A list in the same order as the structure sites containing the summed displacements of each site
    """
    displacement_list = []
    for i in range(len(structure_array[0].sites)):
        displacement_along_path = get_site_displacement(structure_array, i)
        displacement_list.append(displacement_along_path)
    return displacement_list


def get_min_cation_distance(structure, diffuser_index):
    anions = ["O", "S", "F"]
    cation_distances = []
    for i in range(len(structure.sites)):
        if i != diffuser_index:
            if str(structure.sites[i].specie) not in anions:
                cation_distances.append(structure.get_distance(diffuser_index, i))
    return np.amin(cation_distances)


def get_min_anion_distance(structure, diffuser_index):
    anions = ["O", "S", "F"]
    anion_distances = []
    for i in range(len(structure.sites)):
        if i != diffuser_index:
            if str(structure.sites[i].specie) in anions:
                anion_distances.append(structure.get_distance(diffuser_index, i))
    return np.amin(anion_distances)


def get_max_index(list_from_db):
    num_list = []
    for i in range(len(list_from_db)):
        num_list.append(float(list_from_db[i]))
    max_e = max(num_list)
    return num_list.index(max_e)


def find_closest_cation(structures, diff_idx):
    mininmum_distances = []
    for s in range(len(structures)):
        mininmum_distances.append(get_min_cation_distance(structures[s], diff_idx))
    return mininmum_distances


def find_closest_anion(structures, diff_idx):
    mininmum_distances = []
    for s in range(len(structures)):
        mininmum_distances.append(get_min_anion_distance(structures[s], diff_idx))
    return mininmum_distances


def vol_per_anion(structure):
    anions = ["O", "S", "F"]
    num_anions = 0
    for ion in anions:
        num_anions += float(len(structure.indices_from_symbol(ion)))
    return structure.lattice.volume/num_anions


def calculate_average_change(quantity_array):
    average_array = []
    for i in range(len(quantity_array) - 1):
        av_metric = math.sqrt((quantity_array[i] - quantity_array[i + 1])**2)
        average_array.append(av_metric)
    av_metric = np.array(av_metric)
    return np.amax(np.array(quantity_array)) - np.amin(np.array(quantity_array))


def find_intermediate_value(property_array):
    if len(property_array) % 2 == 0:
        intermediate_value = (property_array[len(property_array)/2] + property_array[(len(property_array)/2)-1])/2
    else:
        intermediate_value = property_array[int(math.floor(len(property_array)/2))]
    return intermediate_value


def site_energy_difference(path_energies):
    intermediate_site_energy = find_intermediate_value(path_energies)
    starting_site_energy = path_energies[0]
    return intermediate_site_energy - starting_site_energy


def return_function(property1, property2):
    values = []
    for i in range(len(property1)):
        value = (7.84/property2[i]**2) + (4/property1[i]**2)
        values.append(value)
    return values

# Parse inputs.
# parser = argparse.ArgumentParser()
# parser.add_argument("db_config_path", help="Global path to a YAML file containing user info and database config.")
# args = parser.parse_args()
# yaml_file = str(args.db_config_path)
yaml_file = "/Users/dchannah/yaml/paper2.yaml"
mp_collection = initialize_database(yaml_file)

# User defined:
xproperty = "distance"
yproperty = "path_energy"
search_criterion = " "
search_value = " "
plot_all = False
write_dat = False
color_by_tm = False
spline_flag = False
custom = True
user_query_string = "{\"$or\": [{\"path_id\": \"mp-29356-NEB-001\"}, {\"path_id\": \"mvc-11138-NEB-001\"}, {\"path_id\": \"mvc-12602-NEB-001\"}]}"
# user_query_string = ""

if search_criterion != " " and not custom:
    query_string = "{\"" + search_criterion + "\": \"" + search_value + "\"}"
    query = eval(query_string)
elif custom:
    query = eval(user_query_string)
else:
    query = {}

document_list = mp_collection.find(query)

# Don't mess with this stuff:
x_dictionary = {}
y_dictionary = {}
for item in document_list:
    x_dictionary[item["path_id"]] = get_quantity(item, xproperty)
    y_dictionary[item["path_id"]] = get_quantity(item, yproperty)

document_list = mp_collection.find(query)


if plot_all and color_by_tm:
    plot_collective_tm(x_dictionary, y_dictionary, xproperty, yproperty, document_list)
elif plot_all:
    plot_collective(x_dictionary, y_dictionary, xproperty, yproperty, document_list)
elif color_by_tm:
    plot_individual_tm(x_dictionary, y_dictionary, xproperty, yproperty, document_list, spline_flag)
else:
    plot_individual(x_dictionary, y_dictionary, xproperty, yproperty, document_list, spline_flag)

if write_dat and not plot_all:
    document_list = mp_collection.find(query)
    write_outfile_individual(x_dictionary, y_dictionary, xproperty, yproperty, document_list)

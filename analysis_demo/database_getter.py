import math
import numpy as np
from pymatgen import Structure, Site

"""
This is a collection of methods to get quantities from documents stored in the 
ionic mobility MongoDB.  For somewhat historical reasons, some quantities are 
stored in the database and some are derived on the fly; this may change in the 
future.  These methods could use some refactoring (many repetitive calls to e.g.
get_diffuser_index), but they work for now.  To-do:
    - Add transition metal radius and charge fetching.
    - Other properties...
"""

# Global variables
stored = ["path_energy", "path_length", "coordination", 
          "effective_coordination", "average_bond_length", "distortion_index",
          "activation_energy"]
derived = ["path_percentage", "volume_per_anion", "min_cation_distance", 
           "min_anion_distance", "av_delta_cn", "site_e_diff", 
           "site_cn_diff", "min_cat_dist", "min_anion_dist", 
           "displacements_by_site"]
working_ions = ["Li", "Na", "Ca", "Mg", "Zn"]
anions = ["O", "S", "Se", "F"]


def call_derivation(prop, doc):
    """call_derivation
    A helper method to call individual functions for property derivation.
    :param prop: The property we seek to compute (String)
    :param doc: A MongoDB document.
    :return: Results of calling the specific derivation function.
    """
    func_dict = {
            "path_percentage": path_percentage,
            "volume_per_anion": vol_per_anion,
            "min_cation_distance": min_cation_activated,
            "min_anion_distance": min_anion_activated,
            "av_delta_cn": calc_cn_change,
            "site_e_diff": site_energy_difference,
            "site_cn_diff": site_cn_diff,
            "min_cat_dist": min_cation_along_path,
            "min_anion_dist": min_anion_along_path,
            "displacements_by_site": get_all_displacements,
            }
    return func_dict[prop](doc)


def get_diffuser_index(structures, ion):
    """get_diffuser_index
    From a list of NEB images (Structure objects), finds index of diffuser.
    :param structures: A list of Pymatgen Structure object (NEB trajectory)
    :param ion: Diffusing ion identity (String)
    :return: Integer giving the index of the diffusing atom in species list.
    """
    diff_ion_index, max_path_length = 0.0, 0.0
    for i, s in enumerate(structures[0].sites):
        if str(s.specie) == ion:
            path_length = get_path_length(structures, i)
            if path_length > max_path_length:
                max_path_length = path_length
                diff_ion_index = i
    return diff_ion_index


def get_path_length(structures, d_i):
    """get_path_length
    Gets the path length along the NEB trajectory.
    :param structures: List of structures comprising the NEB trajectory.
    :param d_i: Index of the atom whose diffusion path we seek.
    :return: Path length in angstroms
    """
    d = 0
    for i in range(len(structures) - 1):
        d += structures[i][d_i].distance(structures[i + 1][d_i])
    return d


def determine_working_ion(s):
    """determine_working_ion
    Determines the identity of the working ion.
    :param s: A structure object.
    :return: A String matching the working ion element.
    """
    comp_dict = s.composition.as_dict()
    return [el for el in comp_dict if el in working_ions][0]


def path_percentage(doc):
    """path_percentage
    Normalizes the NEB image coordinates to a distance along the path.
    :param doc: MongoDB document.
    :return: A list of percentage values for each image.
    """
    structures = [Structure.from_dict(s) for s in doc["neb_images"]]
    w_i = determine_working_ion(structures[0])
    d_i = get_diffuser_index(structures, w_i)
    path_total = doc["NEB_analysis"]["path_length"]
    percentage_list = [0]
    for i in range(len(structures) - 1):
        path_distance += structures[i][d_i].distance(structures[i + 1][d_i])
        percentage_list.append(path_distance / path_sum*100)
    return percentage_list


def vol_per_anion(doc):
    """vol_per_anion
    Calculates the volume per anion in a structure.
    :param doc: MongoDB document
    :return: A volume in A^3/anion in unit cell.
    """
    s = Structure.from_dict(doc["neb_images"][0])
    num_anions = 0
    for ion in anions:
        num_anions += len(s.indices_from_symbol(ion))
    return s.lattice.volume/num_anions


def min_cation_activated(doc):
    """min_cation_activated
    Gets the minimum distance to the cation in the activated state.
    :param doc: MongoDB document.
    :return: Minimum activated-state cation distance along path.
    """
    structures = [Structure.from_dict(s) for s in doc["neb_images"]]
    neb_energies = [float(e) for e in doc["NEB_analysis"]["path_energy"]]
    max_energy = max(neb_energies)
    activated_state_index = neb_energies.index(max_energy)
    w_i = determine_working_ion(structures[0])
    d_i = get_diffuser_index(structures, w_i)
    return get_cat_dist_in_struct(structures[activated_state_index], d_i)


def min_anion_activated(doc):
    """min_anion_activated
    Gets the mininmum distance to the anion in the activated state.
    :param doc: MongoDB document.
    :return: Minimum activation state anion distance along the path.
    """
    structures = [Structure.from_dict(s) for s in doc["neb_images"]]
    neb_energies = [float(e) for e in doc["NEB_analysis"]["path_energy"]]
    max_energy = max(neb_energies)
    activated_state_index = neb_energies.index(max_energy)
    w_i = determine_working_ion(structures[0])
    d_i = get_diffuser_index(structures, w_i)
    return get_anion_dist_in_struct(structures[activated_state_index], d_i)


def get_cat_dist_in_struct(s, d_i):
    """get_cat_dist_in_struct
    Finds the nearest distance to a cation in structure s.
    :param s: A Pymatgen Structure object.
    :param d_i: Index of the diffusing ion.
    :return: Distance to the closest cation from the diffusing ion.
    """
    cat_dists = [site.distance(s.sites[d_i]) for i, site in enumerate(s.sites)
                   if i != d_i and str(site.specie) not in anions]
    return min(cat_dists)


def get_anion_dist_in_struct(s, d_i):
    """get_anion_dist_in_struct
    Finds the nearest distance to an anion in structure s.
    :param s: A Pymatgen Structure object.
    :param d_i: Index of the diffusing ion.
    :return: Distance to the closest anion from the diffusing ion.
    """
    anion_dists = [site.distance(s.sites[d_i]) for i, site in enumerate(s.sites)
                   if i != d_i and str(site.specie) in anions]
    return min(anion_dists)


def calc_cn_change(doc):
    """calc_cn_change
    Calculates the change in effective coordination number along the NEB path.
    :param doc: MongoDB document.
    :return: The average change in ECN.
    """
    avgs = []
    ecn_list = np.array(doc["NEB_analysis"]["effective_coordination"])
    max_coord = np.amax(ecn_list)
    min_coord = np.amin(ecn_list)
    return max_coord - min_coord


def site_energy_difference(doc):
    """site_energy_difference
    Returns the energy difference (in meV) between images 0 and N/2.
    :param doc: MongoDB document.
    :return: Out of N images, returns E_{N/2} - E_{0} (in units of meV)
    """
    es = [float(e) for e in doc["NEB_analysis"]["path_energy"]]
    return get_intermediate_difference(es)


def site_cn_diff(doc):
    """site_cn_difference
    Gets the change in coordination number from starting to intermediate state.
    :param doc: MongoDB document.
    :return: ECN_{N/2} - ECN_{0} out of N images.
    """
    ecns = [float(ecn) for ecn in doc["NEB_analysis"]["effective_coordination"]]
    return get_intermediate_difference(ecns)


def get_intermediate_difference(q_a):
    """get_intermediate_difference
    Gets the difference in a quantity between starting and intermediate sites.
    :param q_a: Array of some quantity defined along NEB trajectory.
    :return: Difference in quantity between image 0 and image N/2.
    """
    if len(q_a) % 2 == 0:
        v_int = (q_a[len(q_a)/2] + q_a[(len(q_a)/2) - 1]) / 2
    else:
        v_int = (q_a[int(math.floor(len(q_a)/2))])
    return v_int - q_a[0]


def min_cation_along_path(doc):
    """min_cation_along_path
    Gets the distance to the nearest cation along the whole path.
    :param doc: MongoDB document.
    :return: A list of distances to closest cation along the NEB trajectory.
    """
    structures = [Structure.from_dict(s) for s in doc["neb_images"]]
    neb_energies = [float(e) for e in doc["NEB_analysis"]["path_energy"]]
    max_energy = max(neb_energies)
    activated_state_index = neb_energies.index(max_energy)
    w_i = determine_working_ion(structures[0])
    d_i = get_diffuser_index(structures, w_i)
    cat_dists = [get_cat_dist_in_struct(s, d_i) for i, s in 
                 enumerate(structures)]
    return cat_dists


def min_anion_along_path(doc):
    """min_anion_along_path
    Gets the distance to the nearest anion along the whole path.
    :param doc: MongoDB document.
    :return: A list of distances to closest anion along NEB trajectory.
    """
    structures = [Structure.from_dict(s) for s in doc["neb_images"]]
    neb_energies = [float(e) for e in doc["NEB_analysis"]["path_energy"]]
    max_energy = max(neb_energies)
    activated_state_index = neb_energies.index(max_energy)
    w_i = determine_working_ion(structures[0])
    d_i = get_diffuser_index(structures, w_i)
    anion_dists = [get_anion_dist_in_struct(s, d_i) for i, s in 
                 enumerate(structures)]
    return anion_dists


def get_all_displacements(doc):
    """get_all_displacements
    Get the summed displacement of all sites along the NEB trajectory.
    :param doc: MongoDB document.
    :return: [summed_displacement_site_0, ..., summed_displacement_site_N]
    """
    displacement_list = []
    structures = [Structure.from_dict(s) for s in doc["neb_images"]]
    for i, s in enumerate(structures[0].sites):
        displacement_along_path = get_site_displacements(structures, i)
        displacement_list.append(displacement_along_path)
    return displacement_list


def get_site_displacement(s_a, site_index):
    """get_site_displacement
    Gets the displacement of a particular site along the trajectory.
    :param s_a: Array of structures across the trajectory.
    :param site_index: Index of the site whose displacement we seek.
    :return: Summed displacement along the path.
    """
    # Note: Site.distance returns vector magnitudes, so no worries about signage.
    displacement = 0
    for i in range(len(s_a) - 1):
        site_in_current_img = s_a[i][site_index]
        site_in_next_img = s_a[i + 1][site_index]
        displacement += site_in_current_img.distance(site_in_next_img)
    return displacement

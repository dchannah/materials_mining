#!/usr/bin/python

import numpy as np
import scipy.constants as const
from pymatgen import Composition

# Global vars
EV_PER_ATOM_TO_J_PER_MOL = const.e * const.N_A
ELECTRON_TO_AMPERE_HOURS = EV_PER_ATOM_TO_J_PER_MOL / 3600
allowed_ox = {
        "Ti": [2, 3, 4],
        "V": [2, 3, 4, 5],
        "Cr": [2, 3, 4, 5, 6],
        "Mn": [2, 3, 4, 5],
        "Fe": [2, 3, 4],
        "Co": [2, 3, 4],
        "Ni": [2, 3, 4],
        "Mo": [2, 3, 4, 5, 6]
        }

tms = ["Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Mo"] 
anions = ["O", "S", "Se"]
mobile_ion = ["Li", "Na", "Mg", "Ca", "Zn", "K", "Ag", "Cu"]
mobile_ion_charges = {
        "Li": 1,
        "Na": 1,
        "Mg": 2,
        "Ca": 2,
        "Zn": 2,
        "K": 1,
        "Ag": 1,
        "Cu": 2
        }


# Functions
def compute_tm_ox(comp):
    """compute_tm_ox
    Computes TM oxidation state from a Composition.
    :param comp: Composition dictionary object.
    :return: (float) *Average* oxidation state of TM.
    """
    anion_z = -2
    m_z_list = [mobile_ion_charges[el] for el in comp if el in mobile_ion]
    if len(m_z_list) > 0:
        mobile_z = max(m_z_list)
    else:
        mobile_z = 0
    tm_numbers = [comp[el] for el in comp if el in tms]
    anion_numbers = [comp[el] for el in comp if el in anions]
    mobile_ions = [comp[el] for el in comp if el in mobile_ion]
    total_negative_charge = sum(anion_numbers) * anion_z
    compensating_positive_charge = sum(mobile_ions) * mobile_z
    remaining_charge = total_negative_charge + compensating_positive_charge
    return -1 * remaining_charge/sum(tm_numbers)


def max_tm_ox(comp):
    """max_tm_ox
    Computes the maximum average ox state given a composition.
    :param comp: Composition dictionary object.
    :return: (float) Average max oxidation state.
    """
    tm_max_charges = [max(allowed_ox[el]) for el in comp if el in tms]
    return np.average(tm_max_charges)


def min_tm_ox(comp):
    """min_tm_ox
    Computes the minimum average ox state in a given composition.
    :param comp: Composition dictionary object.
    :return: (float) Average min oxidation state.
    """
    tm_min_charges = [min(allowed_ox[el]) for el in comp if el in tms]
    return np.average(tm_min_charges)


def max_removal(comp):
    """max_removal
    Determines the maximum number of anions removable from comp.
    :param comp: Composition dictionary object.
    :return: (int) The number of working ions we can remove.
    """
    current_tm_ox = compute_tm_ox(comp)
    max_tm = max_tm_ox(comp)
    mobile_z = max([mobile_ion_charges[el] for el in comp if el in mobile_ion])
    delta_z = max_tm - current_tm_ox
    print "TM ox can increase by", delta_z
    # print "Anion ox is", mobile_z
    tm_numbers = sum([comp[el] for el in comp if el in tms])
    mobile_number = sum([comp[el] for el in comp if el in mobile_ion])
    num_remove = (delta_z*tm_numbers) / mobile_z
    return min(num_remove, mobile_number)


def max_insertion(comp, wi):
    """max_insertion
    Determines how many ions we can insert into comp.
    :param comp: A (unintercalated) Composition dictionary object.
    :param wi: A String which is the symbol of the working ion.
    :return: (int) The number of working ions we can insert.
    """
    current_tm_ox = compute_tm_ox(comp)
    min_tm = min_tm_ox(comp)
    mobile_z = mobile_ion_charges[wi]
    delta_z = current_tm_ox - min_tm
    tm_numbers = sum([comp[el] for el in comp if el in tms])
    num_insert = (delta_z*tm_numbers) / mobile_z
    return num_insert


def compute_grav_cap(comp, wi, remove=True):
    """compute_grav_cap
    Computes the gravimetric capacity of a composition.
    :param comp: Composition dictionary object.
    :param wi: String which is the symbol of the working ion.
    :param remove: Are we removing? If not, we're inserting.
    :return: Gravimetric capacity in mAh/g.
    """
    c = comp.as_dict()
    weight = comp.weight
    wi_charge = mobile_ion_charges[wi]
    if remove:
        num_wi = max_removal(c)
    else:
        num_wi = max_insertion(c, wi)
        wi_comp = Composition(wi)
        weight += wi_comp.weight * num_wi
    print num_wi, wi_charge
    cap = num_wi * wi_charge * ELECTRON_TO_AMPERE_HOURS
    grav_cap = cap / (weight/1000)
    return grav_cap

c = Composition("MgMo3P3O13")
wi = "Mg"
print compute_grav_cap(c, wi, remove=True)
# num_insert = max_insertion(c.as_dict(), "Ca")
# print "We can insert", num_insert, "working ions into", c

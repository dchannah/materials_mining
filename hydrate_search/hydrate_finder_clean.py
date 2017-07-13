#!/usr/local/bin/python3
# -*- coding: utf-8 -*-

"""
A program to find compounds which have hydrated phases. Note that because phase
transitions are well known to occur upon hydration, we consider a 'hydrated'
match to be any compound having a composition which is the original formula 
plus an integer number of(H2O) formula units.
"""

import json
from collections import Counter
from pymatgen import MPRester, Composition
from set_generator import generate_query

__author__ = "Daniel Hannah"
__email__ = "dchannah@lbl.gov"

def list_of_keys(comp):
    """
    Returns a List of the keys in the dict representation of a composition.
    """
    return [ele for ele in comp.as_dict()]

def is_hydrate_multiple(original_comp, test_comp):
    """
    Checks if a formula is a hydrated multiple of the input formula.
    :param original_comp: A Pymatgen composition object.
    :param test_comp: A Pymatgen composition
    :return: A boolean.
    """
    for n in range(1, 100):
        if original_comp + n*Composition("H2O") == test_comp:
            return True
    return False


def find_oh_matches(comp):
    """
    Adds O and H to the chemical system and finds the corresponding entries.
    """
    chemsys = [ele for ele in comp.as_dict()]
    chemsys.append('O')
    chemsys.append('H')
    with MPRester() as mpr:
        oh_chemsys = mpr.get_entries_in_chemsys(chemsys)
    hydrate_candidates = [ent for ent in oh_chemsys if "O" in
                          ent.composition.as_dict() and
                          "H" in ent.composition.as_dict()]
    return hydrate_candidates


def hydrated_lists(list_of_entries):
    """
    For a list of entries, finds all the hydrated matches.
    """
    for ent in list_of_entries:
        # ecr = ent.composition.reduced_composition
        ecr = Composition(ent['pretty_formula'])
        oh_entries = find_oh_matches(ecr)
        for oh_ent in oh_entries:
            oh_ecr = oh_ent.composition.reduced_composition
            if is_hydrate_multiple(ecr, oh_ecr):
                print(str(oh_ecr) + " is a hydrate of " + str(ecr))
    return


def check_actual_hydrate(struct):
    """
    Checks if the structure actually contains water molecules.
    """
    for site in struct.sites:
        if str(site.specie) == "O":
           neighbors = struct.get_neighbors(site, 1.5)
           h_sites = [pst for pst in neighbors if str(pst[0].specie) == "H"]
           if len(h_sites) >= 2:
               return True
    return False


def main():
    electrolyte_space = ["Si", "Ge", "P", "Sb", "Zr", "P", "As", "Sn"]
    working_ion = ["Na"]
    anion = ["O"]
    cmpd_list = generate_query(electrolyte_space, working_ion + anion)
    print(len(cmpd_list))
    hydrated_matches = {}
    for cmpd1 in cmpd_list:
        comp1 = Composition(cmpd1['pretty_formula'])
        mpid1 = cmpd1['material_id']
        elist1 = list_of_keys(comp1)
        if mpid1 not in hydrated_matches:
            hydrated_matches[mpid1] = {}
            hydrated_matches[mpid1]['wet_matches'] = []
            hydrated_matches[mpid1]['feature_vector'] = []
            hydrated_matches[mpid1]['formula'] = cmpd1['pretty_formula']
        for cmpd2 in cmpd_list:
            comp2 = Composition(cmpd2['pretty_formula'])
            mpid2 = cmpd2['material_id']
            elist2 = list_of_keys(comp2)
            if is_hydrate_multiple(comp1, comp2):
                with MPRester() as mpr:
                    cmpd2_struct = mpr.get_structure_by_material_id(mpid2)
                if check_actual_hydrate(cmpd2_struct):
                    print(str(mpid2) + " is a real hydrate")
                    wet_dict = {"formula": cmpd2['pretty_formula'], "mpid": mpid2}
                    hydrated_matches[mpid1]['wet_matches'].append(wet_dict)
                else:
                    print(str(mpid2) + " is NOT a real hydrate")
    with open('./ox_without_feature.json', 'w') as f:
        json.dump(hydrated_matches, f) 

if __name__ == "__main__":
    main()
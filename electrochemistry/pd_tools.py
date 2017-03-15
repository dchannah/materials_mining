from pymatgen import MPRester, Composition
from pymatgen.phasediagram.maker import PhaseDiagram
from pymatgen.entries.compatibility import MaterialsProjectCompatibility

"""
Collection of routines to build phase diagrams from the Materials Project.
"""

# Useful global variables.
mp_api = "cVkBCOZrv4TehHfw"  # Replace with your MP API.
dl = ["material_id", "spacegroup", "formation_energy_per_atom",
      "energy_per_atom", "e_above_hull", "pretty_formula"]

def get_min_entry(f):
    """get_min_entry
    Returns the lowest energy entry among entries having same formula.
    :param f: A String containing the chemical composition.
    :return: Entry object for most stable entry out of the set.
    """
    with MPRester(mp_api) as m:
        es = m.get_entries(f, inc_structure=True, property_data=dl)
    if len(es) == 0:
        return None
    else:
        ehull_dict = {}
        for e in es:
            if e.data["e_above_hull"] is not None:
                ehull_dict[e] = e.data["e_above_hull"]
        return min(ehull_dict, key=ehull_dict.get)


def build_corrected_pd(entries):
    """build_corrected_pd
    Builds a PD with entries using Mat.Proj. compatibility.
    :param entries: List of ComputedEntry objects
    :return: A Phase diagram object
    """
    corrected = MaterialsProjectCompatibility().process_entries(entries)
    return PhaseDiagram(corrected)


def get_chemsys(f):
    """get_chemsys
    Returns the entries in a chemical system for formula f.
    :param f: (String) A chemical formula
    :return: A list of entries in the chemical system.
    """
    c = [el for el in Composition(f).as_dict()]
    with MPRester(mp_api) as m:
        c_e = m.get_entries_in_chemsys(c, inc_structure="final", 
                                       property_data=dl)
    return c_e

def get_pruned_chemsys(f):
    """get_pruned_chemsys
    Gets the chemsys but without entries devoid of important data.
    :param f: (String) A chemical formula
    :return: A list of good entries in the chemical system.
    """
    unpruned = get_chemsys(f)
    pruned = [e for e in unpruned if e.data["formation_energy_per_atom"] is not
              None and e.data["e_above_hull"] is not None]
    return pruned

def pd_from_formula(f):
    """pd_from_formula
    Builds a phase diagram containing the compound with formula f.
    :param f: (String) A chemical formula
    :return: A (corrected) PD built around formula f.
    """
    entries = get_chemsys(f)
    return build_corrected_pd(entries)

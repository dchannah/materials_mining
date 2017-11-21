# -*- coding: utf-8 -*-
import pd_tools, sys, os, copy, sympy
from v_calc import build_complete_pd, remove_working_ion

"""
I am currently re-writing the calculator to accont for the simpler set of
voltages spelled out in the latest version of the manuscript. We define here
only 2 distinct voltages:

    1) V_{LDP} = This is the voltage asssuming intercalation into the lowest
    energy discharged polymorph. In general there shold always be a topotactic
    match in this case - we have performed local calculations such that a match
    is always available.

    2) V__{LCP} = This is the voltage resulting from intercalation of the
    ion into the lowest energy charged polymorph. In this case, if no match is
    found OR if a match is found but exhibits an Ehull above 100 meV/atom, a
    hypothetical charged structure is derived with Ehull = 100 meV/atom to
    provide a reasonable upper bound on the intercalation voltage.

"""

__author__ = "Daniel Hannah"
__email__ = "dan@danhannah.site"

# Global variables for chemical systems
MONOVALENTS = ["Li", "Na"]
MULTIVALENTS = ["Ca", "Mg", "Zn"]
WORKING_IONS = MONOVALENTS + MULTIVALENTS
TRANSITION_METALS = ["Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni"]
ANIONS = ["O", "S", "Se"]
CALCDIR = "./conv_inter_calcs"


def find_topo_match(ent, ent_list, hypo=False):
    # Build a structure matcher to use
    s_m = StructureMatcher(primitive_cell=False, scale=True,
                           attempt_supercell=True)

    for e in ent_list:
        if not hypo:
            if sm.fit(ent, e):
                return e
        else:
            if sm.fit(ent, e) and e["e_above_hull"] < CUTOFF:
                return e
            else:
                ent_copy = copy.deepcopy(ent)
                ent_copy["formation_energy_per_atom"] += CUTOFF
                return ent_copy


def generate_cmpds(wis, tms, anions, n_tms):
    for ion in wis:
        for tm in tms:
            for x in anions:
                cmpds.append(ion + tm + str(n_tm) + x + str(n_tm * 2))
    return cmpds


def get_stable_entry(l_e):
    ehull_dict = {}
    for e in l_e:
        if e.data["e_above_hull"] is not None:
            ehull_dict[e] = e.data["e_above_hull"]
    return min(ehull_dict, key=ehull_dict.get)


def get_entries(f_str, pda):
    entries = []
    comp = Composition(f_str).reduced_composition
    for e in pda._pd.all_entries:
        if e.composition.reduced_composition == comp:
            form_e = pda._pd.get_form_energy_per_atom(e)
            e.data["formation_energy_per_atom"] = form_e
            entries.append(e)
    return entries



def main():
    """
    The main method begins by reading in the number of transition metals per
    working ion, which ultimately determines the valency of the reaction. Then
    it goes through a series of function calls to generate the voltage data.
    """

    # Read in number of transition metals per working ion from user.
    n_tm = sys.argv[1]

    # Now get a list of directories containing supplemental local calcs.
    calc_directories = [name for name in os.listdir(CALCDIR) if
                        os.path.isdir(os.path.join(folder, name))]

    # Generate a compound list, using a method for this.
    cmpd_list = generate_cmpds(WORKING_IONS, TRANSITION_METALS, ANIONS, n_tm)

    # For each compound we want to calculate the relevant voltage set.
    for disch_cmpd in cmpd_list:
        chg_cmpd = disch_cmpd[2:]  # Also need corresponding charged compound.
        
        # Get a phase diagram & analyzer including local calculations.
        pd = build_complete_pd(cmpd, calc_directories)
        pda = PDAnalyzer(pd)
        
        disch_entries = get_entries(disch_cmpd)
        ch_entries = get_entries(ch_cmpd)

        """
        First, we'll calculate V_LDP for this compound, which consists of
        finding the lowest energy dischaged product and its topotacic match.
        """
        stable_disch = get_stable_entry(disch_entries)
        ldp_topo_match = find_topo_match(remove_working_ion(stable_disch),
                                            ch_entries)
        ldp_voltages = get_voltage(stable_disch, ldp_topo_match, pda)

        """
        Next, we need to get the LCP voltage. This one is going to be
        trickier, because we need to find the discharged structure that matches
        the stable charged compound while empty.

        This is actually a mess. A Monday project...
        """
        stable_ch = get_stable_entry(ch_entry)
        lcp_topo_empty_match = find_topo_match(stable_ch, [remove_working_ion(e)
                                               for e in disch_entries])
    

        


if __name__ = "__main__":
    main()

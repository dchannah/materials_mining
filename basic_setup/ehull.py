#!/usr/local/bin/python3
# -*- coding: utf-8 -*-

import sys
from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.phasediagram.maker import PhaseDiagram
from pymatgen.phasediagram.analyzer import PDAnalyzer
from pymatgen import MPRester

"""
This script calculates the energy above hull for a system with respect to
phase diagrams built from the Materials Project.
"""

def main():
    """
    Main function.
    """
    # Read the calculation results into a ComputedEntry.
    entl = VaspToComputedEntryDrone().assimilate(sys.argv[1])

    # Check if our local calculation is compatible with MP.
    if MaterialsProjectCompatibility().process_entry(entl) is None:
        print("Calculation not compatible with MP.")
        sys.exit(0)

    # Get other entries sharing a chemical system with the results.
    chemsys = [ele for ele in entl.composition.as_dict()]
    with MPRester() as mpr:
        chemsys_entries = mpr.get_entries_in_chemsys(chemsys)

    # Append our local calculation to the list of entries.
    chemsys_entries.append(entl)

    # Process the entries.
    p_e = MaterialsProjectCompatibility().process_entries(chemsys_entries)

    # Build a phase diagram and an analyzer for it.
    p_d = PhaseDiagram(p_e)
    pda = PDAnalyzer(p_d)

    # Scan stable entries for our calculation.
    for ent in p_d.stable_entries:
        if ent.entry_id is None:
            print(str(ent.composition.reduced_formula) + " 0.0")
            sys.exit(0)

    # Scan unstable entries for our calculation and print decomposition.
    for ent in p_d.unstable_entries:
        if ent.entry_id is None:
            dco, ehull = pda.get_decomp_and_e_above_hull(ent)
            pretty_dc = [("{}:{}".format(k.composition.reduced_formula,
                         k.entry_id), round(v, 2)) for k, v in dco.items()]
            print(str(ent.composition.reduced_formula) + " %.3f" % ehull +
                  " " + str(pretty_dc))

if __name__ == "__main__":
    main()

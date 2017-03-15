#!/usr/bin/python

import os
import sys
from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
from pymatgen.apps.battery.insertion_battery import InsertionElectrode
from pymatgen import MPRester

"""
This script will read a list of topotactically related calculations into an
insertion battery object in Pymatgen.
"""

# Useful global variables
mp_api = "cVkBCOZrv4TehHfw"
wi_ids = {"Li": "mp-135", "Na": "mp-127", "Mg": "mp-153", "Ca": "mp-132", 
          "Zn": "mp-79"}


def read_entry(folder, drone):
    """read_entry
    Reads a VASP entry from a folder.  Must have the XML saved.
    :param folder: Folder containing the VASP calculations results.
    :param drone: A VaspToComputedEntryDrone
    :return: A ComputedStructureEntry object.
    """
    return drone.assimilate(folder)


def build_battery(e_l, wi):
    """build_battery
    Builds an InsertionElectrode from ComputedStructureEntries.
    :param s_l: List of entries for various states of charge.
    :param wi: Charge-carrier ion for the battery. (String)
    :return: An InsertionElectrode object.
    """
    with MPRester(mp_api) as m:
        wi_entry = m.get_entry_by_material_id(wi_ids[wi], inc_structure=True)
    return InsertionElectrode(e_l, wi_entry)


def determine_working_ion(e):
    """determine_working_ion
    Determines the identity of the working ion from a ComputedEntry object.
    :param e: A ComputedEntry object.
    :return: A String representing the working ion.
    """
    wi = None
    for el in e.composition.as_dict():
        if el in wi_ids:
            wi = el
            break
    return wi

if __name__ == "__main__":
    # Build a drone (we do it here to avoid building a new one for every call)
    d = VaspToComputedEntryDrone(inc_structure=True)

    # Read in list of directories containing topotactic states of charge.
    dirlist = sys.argv[1:]
    entry_list = [read_entry(f, d) for f in dirlist]

    # Any entry with more than two ions contains the working ion; find it.
    for entry in entry_list:
        if len(entry.composition.as_dict()) > 2:
            working_ion = determine_working_ion(entry)

    # Can't proceed without a working ion.
    if working_ion is None:
        print("Couldn't find a working ion, quitting.")
        sys.exit(0)

    # Provided we have a working ion, build the battery and print voltages.
    else:
        electrode = build_battery(entry_list, working_ion)
        print(electrode.voltage_pairs)

#!/usr/local/bin/python3
# -*- coding: utf-8 -*-

"""
This is a script to set up relaxations using Materials Project parameters.
"""

import argparse
from pymatgen import Structure
from pymatgen.io.vasp.sets import MPRelaxSet

def setup_vasp(pmg_s, n_e):
    """
    :param pmg_s: Pymatgen Structure object
    :param n_e: Number of electrons to add
    :return: A VaspSet
    """
    nelect = MPRelaxSet(pmg_s).nelect
    if n_e > 0:
        print("###################################")
        print("# BE CAREFUL, YOU ADDED ELECTRONS #")
        print("###################################")
        nelect += n_e
    uis = {
        'NELECT': nelect,
        # 'IVDW': 11
    }

    return MPRelaxSet(pmg_s, user_incar_settings=uis)


def main():
    """
    Main function
    """
    # Process command line arguments
    psr = argparse.ArgumentParser(description="MP relax setup script")
    psr.add_argument('-s', type=str, default="./POSCAR", help='Structure file')
    psr.add_argument('-ne', type=int, default=0, help='# of electrons to add')
    psr.add_argument('-fn', type=str, default=None, help='Folder name for run')
    args = psr.parse_args()

    # Get the structure object
    init_s = Structure.from_file(args.s)

    # Get the foldername
    foldername = "./"
    if args.fn is None:
        foldername += str(init_s.composition).replace(" ", "") + "_relax"
    else:
        foldername += args.fn

    mitset = setup_vasp(init_s, args.ne)
    mitset.write_input(foldername)


if __name__ == "__main__":
    main()

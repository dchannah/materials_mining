#!/usr/local/bin/python3
# -*- coding: utf-8 -*-

"""
This is a script to set up an AIMD simulations using VASP structure file.
"""

import argparse
from pymatgen import Structure
from pymatgen.io.vasp.sets import MITMDSet, MITRelaxSet

def setup_vasp(pmg_s, time_dict, temp_dict, n_e):
    """
    :param pmg_s: Pymatgen Structure object
    :param time_dict: {'timestep': timestep in fs, 'length': total length (fs)}
    :param temp_dict: {'initial': starting temp (K), 'final': final temp (K)}
    :param n_e: Number of electrons to add
    :return: A VaspSet
    """
    nelect = MITRelaxSet(pmg_s).nelect
    if n_e > 0:
        print("###################################")
        print("# BE CAREFUL, YOU ADDED ELECTRONS #")
        print("###################################")
        nelect += n_e
    if time_dict['timestep'] < 2:
        print("dt < 2 fs! Are you sure you wanted a short MD timestep?")
    uis = {
        'NELECT': nelect
    }

    # Need to change default NBLOCK and SMASS behavior if we need temp. scaling
    if temp_dict['initial'] != temp_dict['final']:
        uis['SMASS'] = -1
        uis['NBLOCK'] = 50

    return MITMDSet(pmg_s, temp_dict['initial'], temp_dict['final'],
                    time_dict['length'], time_dict['timestep'],
                    user_incar_settings=uis)


def main():
    """
    Main function
    """
    # Process command line arguments
    psr = argparse.ArgumentParser(description="Script to setup AIMD")
    psr.add_argument('-s', type=str, default='./POSCAR', help='Structure file')
    psr.add_argument('-t0', type=str, required=True, help='Starting T (K)')
    psr.add_argument('-tf', type=str, required=True, help='Final T (K)')
    psr.add_argument('-dt', type=float, default=2, help='Timestep (fs)')
    psr.add_argument('-ne', type=int, default=0, help='# of electrons to add')
    psr.add_argument('-l', type=int, required=True, help='Length of run (fs)')
    psr.add_argument('-fn', type=str, default=None, help='Folder name for run')
    args = psr.parse_args()

    # Get a Pymatgen structure object
    init_s = Structure.from_file(args.s)

    # Generate a default folder name if none is given
    foldername = "./"
    if args.fn is None:
        foldername += str(init_s.composition).replace(" ", "") + "_md"
    else:
        foldername += args.fn

    # Setup temperature ramping info
    temp_data = {'initial': args.t0, 'final': args.tf}

    # Setup timing info
    time_data = {'timestep': args.dt, 'length': args.l}

    # Get a VaspInputSet from these args
    mdset = setup_vasp(init_s, time_data, temp_data, args.ne)

    # Write a VASP run to the specified folder
    mdset.write_input(foldername)


if __name__ == "__main__":
    main()

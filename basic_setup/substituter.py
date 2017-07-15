#!/usr/local/bin/python3
# -*- coding: utf-8 -*-

import json
import os
from argparse import ArgumentParser
from pymatgen import Structure
from pymatgen.transformations.standard_transformations import \
SubstitutionTransformation
from pymatgen.transformations.advanced_transformations import\
EnumerateStructureTransformation

"""
This script carries out structure enumerations.
"""

__author__ = "Daniel Hannah"
__email__ = "dchannah@lbl.gov"


def parse_cmd_line(psr):
    """Handles command-line arguments.

    Args:
        psr (ArgumentParser): Argument parser object to add args to.
    Returns:
        A populated ArgumentParser to access the command line data.
    """
    psr.add_argument('-s', type=str, default="./POSCAR", help='Structure file')
    psr.add_argument('-cell', help='Size of supercell (e.g. 2,2,2)')
    psr.add_argument('-ewald', default=False, help='Rank structures by Ewald?')
    psr.add_argument('-sub', type=json.loads, help='Substitution to perform')
    psr.add_argument('-limit', type=int, default=200, help='Number of structs')
    return psr.parse_args()


def process_arguments(arg):
    """Extracts chemically relevant objects from the command line arguments.

    Args:
        arg: A populated ArgumentParser.  This isn't checked.
    
    Returns:
        A PymatgenStructure object, a supercell array, and whether to use Ewald.
    """
    struct = Structure.from_file(arg.s)
    cell_dim = arg.cell.split(",")
    cell_dim = [int(dim) for dim in cell_dim]
    use_ewald = arg.ewald
    sub_dict = arg.sub
    ext_c = arg.limit
    return struct, cell_dim, use_ewald, sub_dict, ext_c


def apply_substitution(orig_struct, sub):
    """Applies a substitution transformation to the structure.

    Args:
        orig_struct (Structure): A Pymatgen Structure of the original structure.
        sub (dict): {thing_to_sub: {thing_to_sub: new_am1, new_thing: new_amt2}}
    
    Returns:
        A now-substituted structure prior to ordering.
    """
    subber = SubstitutionTransformation(sub)
    return subber.apply_transformation(orig_struct)


def std_enumeration(subbed_struct, cell_dim, coll_limit):
    """Performs an enumeration using the EnumLib library.

    Args:
        subbed_struct (Structure): A structure with some amount of disorder.
        cell_dim (list): A list of maximum a, b, and c cell dimension limits.
        coll_limit (int): Limit to the number of structures produced.
    
    Returns:
        A list of structure dictionaries? I need to test this.
    """
    est = EnumerateStructureTransformation(min_cell_size=1, 
                                           max_cell_size=max(cell_dim))
    if coll_limit == 1:
        return est.apply_transformation(subbed_struct, return_ranked_list=False)
    else:
        return est.apply_transformation(subbed_struct, 
                                        return_ranked_list=coll_limit)
    

def write_structures(list_of_structs):
    """Puts a POSCAR file in a numbered directory for each structure in list."""
    for idx, struct_dic in enumerate(list_of_structs):
        fn = "./" + str(struct_dic['structure'].composition) + "_" + str(idx)
        fn = fn.replace(" ", "")
        os.mkdir(fn)
        struct_dic['structure'].to(fmt='poscar', filename=fn + "/POSCAR")
    return


def main():
    """Main method."""
    parser = ArgumentParser(description="Substitution and ordering script")
    args = parse_cmd_line(parser)
    original, supercell, ewald, substitution, limit = process_arguments(args)
    substituted = apply_substitution(original, substitution)
    struct_results = std_enumeration(substituted, supercell, limit)
    write_structures(struct_results)

if __name__ == "__main__":
    main()
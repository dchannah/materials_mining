#!/usr/bin/python

import sys
import numpy as np
from pymatgen import Molecule, Structure
from pymatgen.analysis.structure_analyzer import OrderParameters


def get_i_in_shell(struct, r_1, r_2):
    """
    Gets all sites within a concentric shell.
    :param struct: Pymatgen Structure object
    :param r_1: Inner radial boundary
    :param r_2: Outer radial boundary
    :return: List of site indices
    """
    sites_in_shell = []
    for i, s in enumerate(struct.sites):
        s_x, s_y, s_z = s.coords[0], s.coords[1], s.coords[2]
        r = np.sqrt(s_x**2 + s_y**2 + s_z**2)
        if r_1 < r <= r_2:
            sites_in_shell.append(i)
    return sites_in_shell


def get_ordering_for_sites(struct, sites):
    """
    Gets the ordering parameter for a list of indices.
    :param struct: Pymatgen Structure object
    :param sites: List of sites for which we want the order parameter.
    :return: Dictionary mapping ordering params to indices.
    """
    bop_types = ["q4"]
    paras = [[]]
    cut_r = 3.0
    # results = {}
    results = []
    bops = OrderParameters(bop_types, paras, cutoff=cut_r)
    for i, s in enumerate(sites):
        op = bops.get_order_parameters(struct, s)
        if op[0] is not None:
            results.append(op[0])
    return np.array(results)


def build_neighbor_dict(struct, r_cut):
    """
    Gets the neighbors within a particular distance of all sites.
    :param struct: Pymatgen structure object
    :param r_cut: Cutoff distance for neighbor searching
    :return: Dictionary mapping indicies to neighbor indices.
    """
    results = {}
    for i, s in enumerate(struct.sites):
        neighbor_list = []
        neighbor_data = struct.get_neighbors(s, r_cut, include_index=True)
        for n in neighbor_data:
            neighbor_list.append(n[2])
        results[i] = neighbor_list
    return results

folder = "/Users/dchannah/work/NU/melt_cdse/traj/split/si_converted/"
# f = sys.argv[1]
# print particle[0]
temps = ["400"]
frames = ["1", "2", "3", "4", "5", "6", "7", "8"]
for t in temps:
    bin_dict = {}
    for u_bs in [10, 15, 20, 25, 35]:
        bin_dict[u_bs] = []
    for f in frames:
        fil = folder + t + "_" + f + ".xyz"
        particle = Molecule.from_file(fil)
        for r_bin in ([0, 10], [10, 15], [15, 20], [20, 25], [25, 35]):
            lower_bound = r_bin[0]
            upper_bound = r_bin[1]
            atom_indices = get_i_in_shell(particle, lower_bound, upper_bound)
            shell = get_ordering_for_sites(particle, atom_indices)
            bin_dict[upper_bound].append(np.average(shell))

for u_bs in [10, 15, 20, 25, 35]:
    print u_bs, str(np.average(bin_dict[u_bs]))
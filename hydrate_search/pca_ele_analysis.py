# -*- coding: utf-8 -*-

"""Principal component analysis of hydrated and non-hydrated structures.

This script analyzes chemical space for areas likely to host undiscovered
phases which can incorporate crystalline water. It does so by utilizing
principal component analysis. Following PCA, an elemental analysis of compounds
inside the hydrated cluster is compared to the distribution of elements in the
whole set.

"""

import json
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pymatgen import Composition
from matplotlib.path import Path
from scipy.spatial import ConvexHull
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

__author__ = "Daniel Hannah"
__email__ = "dchannah@lbl.gov"

def elemental_count(list_of_formulas):
    """Counts the number of compounds containing each element in a list.

    This takes a list of chemical formulas and tallies up the number of
    compounds containing each element. Note that it doesn't concern itself
    with how many atoms of that element are present - only whether the element
    is present or not. Na is excluded because it's in every compound.

    Args:
        list_of_formulas (list): A list of valid chemical formulas.

    Returns:
        A dictionary with elemental symbols as keys and frequencies as values.
    """
    count_dict = {}
    label_dict = {}

    # First we get the absolute counts.
    for formula in list_of_formulas:
        elements_present = list(Composition(formula))
        for ele in elements_present:
            label_dict[ele.number] = str(ele)
            if ele.number not in count_dict and str(ele) != 'Na':
                count_dict[ele.number] = 1
            elif str(ele) != 'Na':
                count_dict[ele.number] += 1

    # Now we scale the values.
    total_elemental_count = sum(count_dict.values())
    for key in count_dict:
        count_dict[key] = count_dict[key]/total_elemental_count

    return count_dict, label_dict

def in_hull(test_point, points, hull):
    """Checks if test_point is in a ConvexHull."""
    hull_path = Path(points[hull.vertices])
    return hull_path.contains_point(test_point)


def read_json_file(filename):
    """Loads a JSON file from a filename.

    Args:
        filname (str): Path to the JSON file.

    Returns:
        A JSON object loaded by the Python json module.
    """
    with open(filename, 'r') as json_file:
        return json.load(json_file)


def check_if_hydratable(dry_dict):
    """Checks if the compound has hydrated formula matches.

    Args:
        dry_dict (dict): A single compound dictionary from the hydrate finder.
    
    Returns:
        True if there are hydrated matches, False otherwise.
    """
    return len(dry_dict['wet_matches']) > 0


def build_pandas_df(coll_dict, fv_list):
    """Reformats dictionary output from the finder for plotting in Pandas.

    This just reformats dictionaries spit out by the hydrate finder code and
    preprocceses them to a format which is easier to handle in pandas/sklearn.

    Args:
        coll_dict (dict): Dictionary in format output by the hydrate finder.
        fv_list (list): A list of descriptors to put into the PCA.

    Returns:
        A dictionary formatted as {mpid: {feature1: value, feature2: value}}
    """
    fv_dict = {}
    for cmpd in coll_dict:
        fv_dict[cmpd] = {}
        for dsc in fv_list:
            fv_dict[cmpd][dsc] = coll_dict[cmpd]['feature_vector'][dsc]
            if check_if_hydratable(coll_dict[cmpd]):
                fv_dict[cmpd]['zclass'] = 'Hydrate-forming'
            else:
                fv_dict[cmpd]['zclass'] = 'Unknown'
    return fv_dict


def main():
    """The main function."""
    # Descriptors to include in principal component analysis.
    descriptor_list = [
        # "mean_atomic_radius",
        "sum_atomic_radius",
        # "mean_electronegativity",
        # "volume_per_atom",
        # "mean_ionic_radius",
        "sum_ionic_radius",
        "sum_row_number",
        # "mean_row_number",
        "sum_group_number"
        # "mean_group_number"
        # "alkali_metal_ratio",
        # "alkali_anion_ratio",
        # "metal_anion_ratio",
        # "coordination_number"
    ]

    json_files = ["./oxide_with_feature_more.json"]

    for j_f in json_files:
        json_data = read_json_file(j_f)
    f_v_set = build_pandas_df(json_data, descriptor_list)

    """
    Simple principal component analysis using Pandas + sklearn
    """
    # Read in a data frame.
    d_f = pd.DataFrame.from_dict(f_v_set, orient='index')
    d_f.reindex_axis(sorted(d_f.columns), axis=1)

    # Set up x and y data.
    x_dat = d_f[descriptor_list].values
    y_dat = d_f['zclass'].values

    # Reduce dimensionality of the data set.
    pca = PCA(n_components=2)
    x_standardized = StandardScaler().fit_transform(x_dat)
    y_transformed = pca.fit_transform(x_standardized)

    # Now we create a dictionary to facilitate elemental analysis after PCA.
    point_ref_dict = {}
    for idx, data_label in enumerate(d_f.index):
        point_ref_dict[str(y_transformed[idx])] = data_label

    # Plot the data in 2-dimensional principal component space.
    with plt.style.context('seaborn-darkgrid'):
        pca1 = 0
        pca2 = 1
        for lab, col in zip(('Unknown', 'Hydrate-forming'), ('green', 'red')):
            plt.scatter(y_transformed[y_dat == lab, pca1],
                        y_transformed[y_dat == lab, pca2], label=lab, c=col)
            if lab == 'Hydrate-forming':
                pts = np.column_stack((y_transformed[y_dat == lab, pca1],
                                       y_transformed[y_dat == lab, pca2]))
                hull = ConvexHull(pts)
                for simplex in hull.simplices:
                    plt.plot(pts[simplex, 0], pts[simplex, 1], 'k-', c=col)

        # Now we get the unknown points within the convex hull.
        all_unk_pts = np.column_stack((y_transformed[y_dat == 'Unknown', pca1],
                                       y_transformed[y_dat == 'Unknown', pca2]))

        unk_pts_in_hull = [u_pt for u_pt in all_unk_pts if
                           in_hull(u_pt, pts, hull)]

        mpids_in_hull = [point_ref_dict[str(p_t)] for p_t in unk_pts_in_hull]

        all_formulas = [json_data[mpd]['formula'] for mpd in d_f.index]
        formulas_in_hull = [json_data[mpd]['formula'] for mpd in mpids_in_hull]

        # Now we generate an elemental histogram for points in the hull & total.
        in_hull_ele_count, labels = elemental_count(formulas_in_hull)
        total_ele_count, dummy_labels = elemental_count(all_formulas)
        total_ele_count_shifted = {}
        shift = 0.1
        for key in total_ele_count:
            total_ele_count_shifted[key + shift] = total_ele_count[key]


        # Show the PCA plot
        plt.xlabel("Principal Component 1")
        plt.ylabel("Principal Component 2")
        plt.legend(frameon=True)
        plt.show()

        # Now we look at the histogram
        plt.gcf().clear()  # Clear the old plot
        plt.bar(list(in_hull_ele_count.keys()), in_hull_ele_count.values(),
                color='r', width=0.1, label='Inside hydrate cluster')
        plt.bar(list(total_ele_count_shifted.keys()),
                total_ele_count_shifted.values(), color='g', width=0.1,
                label='All compounds')
        bar_plot_labels = [labels[value] for value in in_hull_ele_count.keys()]
        plt.xticks(list(in_hull_ele_count.keys()), bar_plot_labels)
        plt.ylabel("Fraction of material class (%)")
        plt.legend(frameon=True)
        plt.show()

if __name__ == "__main__":
    main()

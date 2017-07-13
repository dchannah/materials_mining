# -*- coding: utf-8 -*-

"""F-testing tool for descriptors.

This script builds on the existing data-mining tools for hydrate finding
and carries out f_testing of feature vectors.

"""

import json
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

__author__ = "Daniel Hannah"
__email__ = "dchannah@lbl.gov"


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


def split_feature_lists(c_dict, feature):
    """Creates lists of hydrated and non-hydrated feature values for F-testing.

    Note:
        This method sloppily assumes the feature vector field has been filled.
    
    Args:
        c_dict (dict): A dictionary (with dry cmpd keys) of hydrate finds.
        feature (string): Feature whose values we want to split up.

    Returns:
        A dictionary of two lists containing dry and wet values.
    """
    value_dict = {"dry": [], "wet": []}
    for d_c in c_dict:
        f_v_value = c_dict[d_c]['feature_vector'][feature]
        
        # Add the value to the appropriate list
        if check_if_hydratable(c_dict[d_c]):
            value_dict['wet'].append(f_v_value)
        else:
            value_dict['dry'].append(f_v_value)
    return value_dict


def f_value(group1, group2):
    """Carries out an F-test (or similar) using two lists.

    This method compares the variance within hydrates/dry compounds to the
    variance across those groups in order to assess a descriptor's utility.

    Args:
        group1 (list): List of values for the first group (ex. hydrates)
        group2 (list): List of values for the second group (ex. dry compounds)

    Returns:
        The F-ratio or Kruskal-Wallis ratio value.
    """
    return stats.kruskal(group1, group2)
    # return stats.f_oneway(group1, group2)


def main():
    """The main function."""
    descriptor_list = [
        "mean_atomic_radius",
        "sum_atomic_radius",
        "mean_electronegativity",
        "volume_per_atom",
        "mean_ionic_radius",
        "sum_ionic_radius",
        "sum_row_number",
        "mean_row_number",
        "sum_group_number",
        "mean_group_number",
        "coordination_number"
    ]
    json_files = ["oxide_with_feature_more.json"]

    for j_f in json_files:
        json_data = read_json_file(j_f)
        bar_plot_labels = []
        f_values = []
        for descriptor in descriptor_list:
            bar_plot_labels.append(descriptor)
            data_dict = split_feature_lists(json_data, descriptor)
            f_values.append(f_value(data_dict['dry'], data_dict['wet'])[0])

    bar_y_positions = np.arange(len(f_values))
    plt.bar(bar_y_positions, f_values)
    plt.xticks(bar_y_positions, bar_plot_labels, rotation=90)
    plt.ylabel("Kruskal-Wallis")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()

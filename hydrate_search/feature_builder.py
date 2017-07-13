# -*- coding: utf-8 -*-

"""The feature builder.

This program populates the "feature_vector" field of the dictionary objects
spit out by our hydrate-finding code and proceeds to do some learning and
clustering analysis on them.

"""

import json
from new_feature_finder import FeatureFinder

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


def populate_feature_vectors(f_v, c_dict):
    """Populates all of the feature vector fields in a compound collection.

    This method takes a dictionary of dictionaries relating dry compounds to
    their hydrated matches (if any). It creates a feature vector based on the
    requested features (a list passed to to the method) and stored in the
    "feature_vector" field of the dictionary.

    Args:
        f_v (list) = A list of descriptors to include in the feature vector.
        c_dict (dict) = A dictionary (with dry cmpd keys) of hydrate finds.

    """
    for d_c in c_dict:
        c_dict[d_c]["feature_vector"] = populate_f_v(f_v, d_c, c_dict[d_c])
    return


def populate_f_v(vector, dry_id, dry_dict):
    """Creates a feature vector for a dry compound.

    This method creates a FeatureFinder object and uses it to generate the
    feature vector for the dry compound in question based on a user-supplied
    list of descriptors.

    Args:
        vector (list): A List of descriptors to include in the vector.
        dry_id (str): Materials Project ID of the dry compound.

    Returns:
        A dictionary of descriptors to build the feature vector.
    """
    f_f = FeatureFinder(dry_id, dry_dict)
    f_v = {}
    for descriptor in vector:
        f_v[descriptor] = f_f.feature_finder(descriptor)
    return f_v


def main():
    """The main function."""
    # We define the descriptor list here to make this easy to read from a file.
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
        "alkali_metal_ratio",
        "alkali_anion_ratio",
        "metal_anion_ratio",
        "coordination_number"
    ]

    # A list of JSON files from our hydrate finder.
    json_files = ["sample_json_no_feature.json"]

    for j_f in json_files:
        json_data = read_json_file(j_f)
        populate_feature_vectors(descriptor_list, json_data)
        with open('./sample_json_with_feature.json', 'w') as f:
            json.dump(json_data, f)

if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-

"""Frequency examination for variables.

This is to check whether the features we've chosen for our analysis are
distributed normally across the two groups (hydrated and dry compounds).

"""

from f_tester import split_feature_lists, read_json_file
import matplotlib.pyplot as plt


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
        dict_of_data_dict = {}
        for descriptor in descriptor_list:
            bar_plot_labels.append(descriptor)
            data_dict = split_feature_lists(json_data, descriptor)
            dict_of_data_dict[descriptor] = data_dict
    
    # Visualize histograms of the data
    fig, axarr = plt.subplots(2, len(dict_of_data_dict))
    for idx, d_e in enumerate(dict_of_data_dict):
        # Plot the histograms
        axarr[0, idx].hist(dict_of_data_dict[d_e]['dry'], color='red')
        axarr[1, idx].hist(dict_of_data_dict[d_e]['wet'], color='blue')

        # Hide the y-axes (we are just checking the distribution)
        axarr[0, idx].axes.get_yaxis().set_visible(False)
        axarr[1, idx].axes.get_yaxis().set_visible(False)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.17)  # Make readable on laptop

    plt.show()

if __name__ == "__main__":
    main()
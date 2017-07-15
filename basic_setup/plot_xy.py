#!/usr/local/bin/python3
# -*- coding: utf-8 -*-

"""
This is a general purpose x vs. y plotting script with decent-looking defaults.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.texmanager

"""
Global user-defined variables below; user should change only these.
"""
WANT_SCATTER = True  # Do you want a scatter plot?
WANT_LINES = False  # Do you want a line plot (w/ linear interpolation?)
STYLE_DICT = {  # Style settings.
    'lw': 2.5,  # Linewidth
    # 'labels': None,  # If used: [file1_label, file2_label, ..., fileN_label]
    'labels': ["2", "3", "4 (tet)", "4 (sq)", "5", "6", "8"],
    'size': 85,  # Marker size for color map.
    'ec': 'black'  # Edge color of scatter plot markers.
}
COLUMNS = {'x': 0, 'y': 1}  # Column number locations
CMAP = plt.cm.jet
ANNOTATIONS = []  # Text labels (in order) for data points
DRAW_HORIZONTAL_LINES = False
DRAW_VERTICAL_LINES = False
HORIZONTAL_LINES = [0]
VERTICAL_LINES = [0.0313]
HIDEX = False  # Hide x-axis?
HIDEY = False  # Hide y-axis?

"""
Stuff that users shouldn't need to change.
"""

def add_scatter_plot(axi, xdat, ydat, style_params):
    """
    Adds a scatter plot to an axis object.
    :param axi: A matplotlib Axis.
    :param xdat: Data for the x-axis.
    :param ydat: Data for the y-axis.
    :param style_params: Dict of lw, color, size, label
    :return: None
    """
    axi.scatter(xdat, ydat,
                zorder=2,
                linewidth=style_params['lw'],
                edgecolor=style_params['ec'],
                color=style_params['color'],
                s=style_params['size'])
    if len(ANNOTATIONS) == len(xdat):
        for i, txt in enumerate(ANNOTATIONS):
            axi.annotate(txt, (xdat[i], ydat[i]), xytext=(-10, 10),
                         textcoords='offset points')
    return


def add_line_plot(axi, xdat, ydat, style_params):
    """
    Adds a line plot to an axis object.
    :param axi: A matplotlib Axis.
    :param xdat: Data for the x-axis.
    :param ydat: Data for the y-axis.
    :param style_params: Dict of lw, color, size, label
    :return: None
    """
    axi.plot(xdat, ydat,
             linewidth=style_params['lw'],
             zorder=1,
             color=style_params['color'],
             label=style_params['label'])
    return


def plot_files(axi, file_list, column_nos, plot_types, style_params):
    """
    Reads files and adds them to the plot.
    :param axi: Matplotlib axis object
    :param file_list List of files we seek to plot.
    :param column_nos: {'x': Col. no for x data, 'y': Col no. for y data}
    :param plot_types: {'scatter': True/False, 'lines': True/False}
    :param style_params: Style dictionary (lw, color, size, label, etc.)
    :return: None
    """
    # Setup color map
    cmaplist = [CMAP(i) for i in range(CMAP.N)]  # Color map

    # Setup even sampling of color map based on number of files.
    spacing = int(len(cmaplist) / len(file_list))

    for i, datafile in enumerate(file_list):
        style_params['color'] = cmaplist[i * spacing]
        if style_params['labels'] is not None:
            style_params['label'] = style_params['labels'][i]
        else:
            style_params['label'] = ""
        print("Trying to read " + str(datafile))
        with open(datafile, 'r') as file:
            data = np.loadtxt(file)
        x_data = data[:, column_nos['x']]
        y_data = data[:, column_nos['y']]
        if plot_types['scatter']:
            add_scatter_plot(axi, x_data, y_data, style_params)
        if plot_types['lines']:
            add_line_plot(axi, x_data, y_data, style_params)
    return


def draw_lines(axi, which_lines, line_dict):
    """
    Draws dashed lines on the plot.
    :param axi: Matplotlib axis object
    :param which_lines: {'horiz': True/False, 'vert': True/False}
    :param line_dict: {'horiz': [line1, ..., lineN], 'vert': [line1, ..., N]}
    :return: None
    """
    lstyle = 'dashed'
    clr = 'black'
    if which_lines['horiz']:
        for y_loc in line_dict['horiz']:
            axi.axhline(y=y_loc, linestyle=lstyle, color=clr)
    if which_lines['vert']:
        for x_loc in line_dict['vert']:
            axi.axvline(x=x_loc, linestyle=lstyle, color=clr)
    return


def make_pretty_plot(xlabel, ylabel):
    """
    Makes a nice looking plot to work on.
    :param xlabel: x-axis label for plot
    :param ylabel: y-axis label for plot
    :return: Matplotlib axis object.
    """
    # Create the figure and axis.
    fig = plt.figure(figsize=(6, 6))
    axi = fig.add_subplot(111)

    # Hide axes if needed.
    if HIDEX:
        axi.axes.get_xaxis().set_visible(False)
    if HIDEY:
        axi.axes.get_yaxis().set_visible(False)

    # Thick borders
    border_width = 3
    for axis in ['top', 'bottom', 'left', 'right']:
        axi.spines[axis].set_linewidth(border_width)

    # Label and font sizes
    axi.set_xlabel(xlabel, fontsize=24)
    axi.set_ylabel(ylabel, fontsize=24)
    axi.tick_params(axis='x', labelsize=22)
    axi.tick_params(axis='y', labelsize=22)
    return axi


def main():
    """
    Reads the files and axis labels in from argparse and then plots them.
    """
    # Read in the arguments from the command line.
    prs = argparse.ArgumentParser(description='General plotting script.')
    prs.add_argument('-fl', nargs='+', help='List of files to plot.')
    prs.add_argument('-x', type=str, default='x', help='x label')
    prs.add_argument('-y', type=str, default='y', help='y label')
    args = prs.parse_args()

    # Set up a square figure using matplotlib.
    ax1 = make_pretty_plot(args.x, args.y)

    # Plot files on this axis.
    what_to_plot = {'scatter': WANT_SCATTER, 'lines': WANT_LINES}
    plot_files(ax1, args.fl, COLUMNS, what_to_plot, STYLE_DICT)

    # Draw any needed lines on the plot
    draw_lines(ax1,
               {'horiz': DRAW_HORIZONTAL_LINES, 'vert': DRAW_VERTICAL_LINES},
               {'horiz': HORIZONTAL_LINES, 'vert': VERTICAL_LINES})

    # Tighten layout, add legend if needed, and show plot.
    if STYLE_DICT['labels'] is not None:
        plt.legend(loc='best')
    plt.tight_layout()
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()

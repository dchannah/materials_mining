#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt


def autolabel(sp, rects, top, bottom, names):
    """autolabel
    Adds specified text labels to bars.
    :param sp: The particular subplot we want to label
    :param rects: A barh object
    :param top: The maximum voltage - for a horizontal plot, 'right'
    :param bottom: The minimum voltage - for a horizontal plot, 'left'
    :param names: The array of labels
    """
    tol = 0.085
    for ii,rect in enumerate(rects):
        if np.abs(top[ii] - bottom[ii]) < tol:
            higherval = max(top[ii], bottom[ii])
            xloc = higherval + 0.02
            halign = 'left'
        else:
            xloc = (top[ii] + bottom[ii])/2
            halign= 'center'
        yloc = rect.get_y() + rect.get_height()/2
        sp.text(xloc, yloc, names[ii], horizontalalignment=halign, verticalalignment='center', fontsize=16)
    return


def autolabel_side(sp, rects, top, bottom, names):
    """autolabel_side
    Generates side labels for the bars.  
    :param sp: The particular subplot we want to label
    :param rects: A barh object
    :param top: The maximum voltage - for a horizontal plot, 'right'
    :param bottom: The minimum voltage - for a horizontal plot, 'left'
    :param names: The array of labels
    """
    for ii, rect in enumerate(rects):
        if "Zn" not in names[ii]:
            baseval = min(top[ii], bottom[ii])
            xshift = -0.02
            xloc = baseval + xshift
            halign = 'right'
            yloc = rect.get_y() + rect.get_height()/2
            sp.text(xloc, yloc, names[ii], horizontalalignment=halign, verticalalignment='center', fontsize=16)
        else:
            baseval = max(top[ii], bottom[ii])
            xshift = 0.02
            xloc = baseval + xshift
            halign = 'left'
            yloc = rect.get_y() + rect.get_height()/2
            sp.text(xloc, yloc, names[ii], horizontalalignment=halign, verticalalignment='center', fontsize=16)
    return

def generate_colors(delta_v):
    """generate_colors
    Adds a list of color instructions based on preference for intercalation over conversion.
    :param delta_v: v_int - v_conv
    """
    colors = []
    for i, dv in enumerate(delta_v):
        if dv < 0:
            colors.append('green')
        else:
            colors.append('red')
    return colors

# Need a dictionary of working ions
working_ions = ["Li", "Na", "Mg", "Ca", "Zn"]
results_dict = {}
for i in working_ions:
    results_dict[i] = []

# Create stacked plots to plot by working ion
fig, axarr = plt.subplots(len(working_ions), 1, sharex=True)
subplots_dict = {}
for i, wi in enumerate(working_ions):
    subplots_dict[wi] = axarr[i]

# Read in the data as labels and voltages
with open(sys.argv[1], 'r') as f:
    label_data = np.genfromtxt(f, usecols=[0], dtype='string')
with open(sys.argv[1], 'r') as f:
    num_data = np.genfromtxt(f, usecols=[1, 2, 3, 4])

# Partition the data up into the results dictionary
for i, l in enumerate(label_data):
    wi_label = l[0:2]
    if wi_label not in working_ions:
        print "Compound name is messed up, not counting", l
    else:
        results_dict[wi_label].append([l, num_data[i]])

# Plot by working ion
for wi in working_ions:
    subplot = subplots_dict[wi]
    plt_labels = []
    plt_set_max = []
    plt_set_min = []
    plt_set_max_lefts = []
    plt_set_min_lefts = []
    c_max_conv = []
    c_min_conv = []
    for c in results_dict[wi]:
        label = c[0]
        v = c[1]
        min_int = v[0]
        max_int = v[1]
        min_conv = v[2]
        max_conv = v[3]
        max_window = max_int - min_int
        min_window = max_conv - min_conv
        max_left = min(max_int, min_int)
        min_left = min(max_conv, min_conv)
        plt_set_max.append(max_window)
        plt_set_min.append(min_window)
        plt_set_max_lefts.append(max_left)
        plt_set_min_lefts.append(min_left)
        plt_labels.append(label)
        c_max_conv.append(max_int)
        c_min_conv.append(min_int)
    vbar_centers = np.arange(len(plt_labels))
    # color_max = generate_colors(plt_set_max)
    # color_min = generate_colors(plt_set_min)
    color_max, color_min = "green", "red"
    rect_max = subplot.barh(vbar_centers, np.absolute(plt_set_max), left=plt_set_max_lefts, alpha=0.5, color=color_max, label="Intercalation")
    rect_min = subplot.barh(vbar_centers, np.absolute(plt_set_min), left=plt_set_min_lefts, alpha=0.25, hatch='//', color=color_min, label="Conversion")
    subplot.text(0.075, 0.85, wi, horizontalalignment='center', verticalalignment='center', transform=subplot.transAxes, fontsize=24)
    autolabel_side(subplot, rect_max, c_max_conv, c_min_conv, plt_labels)

# Label only the bottom axis
subplots_dict["Zn"].set_xlabel("Voltage (V)", fontsize=26)
subplots_dict["Zn"].tick_params(axis='x', labelsize=24)
subplots_dict["Zn"].legend(loc='best')
for ax in subplots_dict:
    subplots_dict[ax].axes.get_yaxis().set_visible(False)
    subplots_dict[ax].grid()


# Making it fancy
plt.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',
       r'\sisetup{detect-all}',
       r'\usepackage{helvet}',
       r'\usepackage[eulergreek,EULERGREEK]{sansmath}',
       r'\sansmath'
]
border_width = 2
for ax in subplots_dict:
    [i.set_linewidth(border_width) for i in subplots_dict[ax].spines.itervalues()]
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')

# Plot the thing!
plt.tight_layout
plt.show()

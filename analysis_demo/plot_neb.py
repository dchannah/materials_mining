#!/usr/bin/python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

# Read all files into a list for plotting.
files = []
for i in range(1, len(sys.argv)):
    with open(sys.argv[i], 'r') as f:
        files.append(np.genfromtxt(f, skip_footer=1))

# Figure settings
fig = plt.figure(figsize=[6,6], dpi=100)
ax = fig.add_subplot(111)


# For multi-column files
x_column_no = 0
y_column_no = 1

# Dictionary objects for labeling and coloring
legend_dict = {0: "Charged", 1: "Discharged"}
color_dict = {0: 'blue', 1: 'blue', 2: 'green', 3: 'orange'}
linestyle_dict = {0: 'solid', 1: 'dashed'}

# Read in and plot all files
for i in range(0, len(files)):
    ax.scatter(files[i][:,x_column_no], files[i][:,y_column_no], s=70, linewidth=2.5, edgecolor='black', zorder=2,
               color=color_dict[i]) 
    # label_plot = ("\\begin{verbatim}" + str(sys.argv[i+1]) + "\\end{verbatim}")
    label_plot = legend_dict[i]
    label_plot = ""
    
    # Uncomment below for spline fit:
    tck = interpolate.splrep(files[i][:, x_column_no], files[i][:, y_column_no], s=0)
    xnew = np.arange(0, 100, 0.1)
    splfit = interpolate.splev(xnew, tck, der=0)
    ax.plot(xnew, splfit, linewidth=2.0, zorder=1, color=color_dict[i], linestyle=linestyle_dict[i], label=label_plot)

    # Uncomment below for standard line plot:
    # ax.plot(files[i][:, x_column_no], files[i][:,y_column_no], linewidth=2, zorder=1, color=color_dict[i], label=label_plot)
# Label axes
label_dict = {0: "Distance along path (\%)", 1: "Relative Energy (meV)", 2: "Bond Valence", 3: "Coordination Number",
              4: "Bond Length", 5: "Effective Coord.", 6: "Distortion Index"}
xlabel = "\\textbf{" + label_dict[x_column_no] + "}"
ylabel = "\\textbf{" + label_dict[y_column_no] + "}"

# Make plot look pretty
plt.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]
ax.set_xlabel(xlabel, fontsize=24)
ax.set_xlim([0, 100])
ax.set_ylabel(ylabel, fontsize=24)
ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)
border_width = 2
[i.set_linewidth(border_width) for i in ax.spines.itervalues()]
plt.tight_layout()
plt.legend(loc='best', prop={'size':12})
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
plt.tight_layout()

# If you want to save a figure:
# output_filename = "MgTiS2.pdf"
# plt.savefig(output_filename)

# Pop-up plot window
plt.show()
# plt.savefig("figure.png")
# Save a PNG file
# filename = os.path.splitext(sys.argv[1])[0]
# plt.savefig(filename + ".svg")

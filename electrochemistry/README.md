# electrochemistry

These tools are for carrying out electrochemical analyses (voltage, mostly) of VASP calculations as well as materials data stored in the Materials Project database.  Scripts:

* battery_builder: A simple script for generating voltage curves and other data (e.g. capacity, volume change, etc.) given DFT calculations of various charge states for an intercalation battery cathode material.

* v_calc.py: Analyses the range of conversion and intercalation voltages possible for all ABX<sub>2</sub> and AB<sub>2</sub>X<sub>4</sub> polymorphs in the Materials Project database, where A is a common battery working ion, B is a 3d transition metal, and X is an anion.  Reads in supplementary local calculations when needed.

* pd_tools.py: A collection of fancier methods which basically call MPRester to do useful things for large scale stability and electrochemistry analyses.

* oxi_filter.py: A simple script to ballpark oxidation states in cathode materials based on "common" oxidation states for each element. Assumes all redox chemistry occurs on the transition metal site. 

* parse\_rxn\_table.py: Converts the output of v\_calc.py into LaTeX tables for easy typesetting (one\_electron_reactions.txt is a sample input file).

* formula\_to\_tex.py: Heavy-lifting for the text parsing in parse\_rxn\_table.py. 
# analysis_demo

This an analysis demo for the ionic mobility database (these are slightly modified version of the tools in my [ionic_mobility](https://github.com/dchannah/materials_mining/tree/master/ionic_mobility) folder).  Two scripts are present here:

* database_getter.py: A collection of methods to pull various stored and derived properties from a MongoDB.

* analyze_db_visualize.py: A script which uses database_getter to query a particular list of compounds from the MongoDB (in the example, we query for all Ca-ion battery cathodes) and an associated set of properties, then regresses the activation energy against each of these properties and displays the errors.

* plot_neb.py: A dynamic, messy script that I use to plot results of individual ionic mobility calculations.  Uses matplotlib and LaTeX for labels.

* neb_profile.dat: The output from NEB analysis; this is fed to plot_neb.py to plot (for example) distance vs. energy.  Example:  
```python plot_neb.py neb_profile.dat ```

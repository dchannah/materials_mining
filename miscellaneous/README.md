# Miscellaneous
Some misc. scripts (and sample data) to facilitate my research.
* plot_results_bar.py: Reads in voltage window data from an external file (sample data in 1e_oxides_bar.dat) and plots conversion and intercalation voltage windows.  1e_oxides_barplot.dat is a sample data set for this script.

* ml_testing.py: A simple machine learning demo to correlate band gaps with chemical composition.  Uses the Materials Project database.

* set_generator.py: Generates Mongo queries to filter materials for further high-throughput computational screening.

* murnaghan.py: Generates bulk moduli by fitting volume-depedendent energy points to the [Murnaghan equation of state](https://en.wikipedia.org/wiki/Murnaghan_equation_of_state).  Uses numpy's polyfit feature to generate initial guesses for fitting.

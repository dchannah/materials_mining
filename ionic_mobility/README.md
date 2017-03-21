# ionic_mobility

This collection of tools is designed to interact with a MongoDB containing the results of Nudged Elastic Band calculations.  Note that these Nudged Elastic Band calculations can use energies derived from anything (electrostatics, density functional theory, etc.) so long as image-wise energies are resolved and stores in the database.

* database_getter.py: Collection of methods to pull various quantities from the ionic mobility database.

* analyze_db.py: Methods to facilitate the generation of feature vectors from the ionic mobility database; also contains other methods useful for querying the DB.

* crm.py: An attempt to apply clustering-ranking-modeling to ionic mobility data.  Currently a work in progress; right now I dump the data and use regress_feature to get the cross_Validation scores for linear regression manually (because there are a small number of features), but the next step is writing a routine to regress all the features from the feature vector in a given cluster against the ionic migration barrier and select the best one.

* database_scanner_clean.py: An older, messier version of database_getter. When I began to explore more data mining approaches I took a more modular approach and decided to rebuild this script in the form of database getter.  This script is being phased out of my workflow but still contains some unique functionality.

* diffusion_updater_v3.py: This script parses data from ionic mobility simulations (namely, Nudged Elastic Band calculations run in VASP) and processes it a bit to get it in shape for the ionic mobility MongoDB.  Contains various useful routines for calculating physical properties (like bond distance) along the NEB trajectory, as well as routines for visualization.

* neb_push_v3.py: Uses the UpdaterNEB class (see diffusion_updater_v3.py) to process and merge NEB trajectories into the ionic mobility MongoDB.

* pca_test.py: Another work in progress - attempting to apply principal component analysis to feature vectors generated from NEB trajectories.   

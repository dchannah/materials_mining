import itertools
import numpy as np
from pymatgen import MPRester, Composition, Element
from pymatgen.core import periodic_table
from sklearn import linear_model, cross_validation, metrics, ensemble

"""
Here, we're just messing around with simple machine learning tools from MP.
"""

mp_api = "cVkBCOZrv4TehHfw"
dl = ["pretty_formula", "band_gap", "e_above_hull"]

MAX_VECTOR = 100


def build_physical_features(comp_set):
    """build_physical_features
    Builds a set of physically motivated features from a set of compositions.
    :param comp_set: A list of strings of chemical compositions.
    :return: A list of lists of physical features.
    """
    vector_list = []
    for comp in comp_set:
        feature_list = []
        atomic_no = []
        eneg = []
        fraction = []
        group = []
        for el in comp:
            atomic_no.append(float(el.Z))
            eneg.append(el.X)
            fraction.append(comp.get_atomic_fraction(el))
            group.append(float(el.group))
        reverse = (fraction[1] > fraction[0])
        for features in [atomic_no, eneg, fraction, group]:
            if reverse:
                features.reverse()
        feature_list.append(eneg[0] - eneg[1])
        feature_list.append(fraction[0] - fraction[1])
        feature_list.append(group[0])
        feature_list.append(group[1])
        vector_list.append(feature_list)
    return vector_list




def get_formula_and_bandgap(chemsys):
    """get_formula_and_bandgap
    Gets the formula and band gap from Materials Project given a chemsys.
    :param chemsys: (String) A binary compound
    :return: [pretty_formula, band_gap]
    """
    binary_cmpd = chemsys[0] + chemsys[1]
    stable_results = []
    with MPRester(mp_api) as m:
	ents = m.get_entries(binary_cmpd, property_data=dl)
    for e in ents:
	if e.data['e_above_hull'] < 1e-6:  # Want only thermodynamically stable.
	    stable_results.append([e.data['pretty_formula'], e.data['band_gap']])
    return stable_results


def naive_vectorize(comp):
    """naive_vectorize
    Converts a chemical composition into numbers based on atomic fractions.
    :param comp: Composition
    :return: A vector based on the composition's stoichiometry.
    """
    comp = Composition(comp)
    vector = np.zeros(MAX_VECTOR)
    for el in comp:
        frac = comp.get_atomic_fraction(el)
	vector[el.Z - 1] = frac  # Python indexing starts at 0, science doesn't
    return vector


def all_symbols(max_nuc=100):
    """all_symbols
    Gets the atomic symbol from the nuclear charge Z.
    :param max_nuc: Biggest nucleus to get.  Defaults to the whole PT.
    :return: A tuple of all the requested symbols corresponding to Zs.
    """
    symbols = []
    for z in range(1, max_nuc + 1):
        s = periodic_table.Element.from_Z(z)
        if s is None:
            break  # Can't have a None in our string list.
        symbols.append(str(s))

    return tuple(symbols)

# First, generate a list of all possible binary compounds in the periodic table.
all_binaries = itertools.combinations(all_symbols(), 2)

# Now we create two lists: Chemical compounds and band gaps.
materials = []
bandgaps = []
naive_features = []

for binary in all_binaries:
    queried_data = get_formula_and_bandgap(binary)
    for stable_compound in queried_data:
        materials.append(stable_compound[0])
	bandgaps.append(stable_compound[1])
        naive_features.append(naive_vectorize(stable_compound[0]))

# If our algorithm doesn't produce better predictions than taking the average,
# we're in trouble.
baseline_e = mean(abs(mean(bandgaps)) - bandgaps)

linear = linear_model.Ridge(alpha = 0.5)  # Linear least squares.
cv = cross_validation.ShuffleSplit(len(bandgaps), n_inter=10, test_size=0.1, 
                                   random_state=0)  # Training & test data
scores = cross_validation.cross_val_score(linear, naiveFeatures, bandgaps,
                                          cv=cv, scoring='mean_absolute_error')

print("Mean absolute error of linear ridge regression on band gap model:")
print(str(round(abs(mean(scores)), 3)) + " eV")
print("For comparison, mean absolute error of just guessing avg. band gap is:")
print(str(round(baseline_e, 3)) + " eV")

# Some more detailed info about the model:
print("Coefficients for each element:")
linear.fit(naive_features, bandgaps)
for i in range(MAX_VECTOR):
    element = Element.from_Z(i + 1)
    print(element.symbol + ": " + str(linear.coef_[i]))

# Instead of just mining the compositions, can we do better with an informed set
# of features?
phys_features = build_physical_features(materials)
phys_scores = cross_validation.cross_val_score(linear, phys_features, bandgaps,
                                               cv=cv, 
                                               scoring='mean_absolute_error')

"""
There is no (physical/chemical) reason to believe band gap should be linearly 
correlated with either composition or our physical features.  So, we will try
a nonlinear regression model - random forest regression.  
"""

rfr = ensemble.RandomForestRegressor(n_estimators=10)
rfr_naive_scores = cross_validation.cross_val_score(rfr, naive_features, 
                                                    bandgaps, cv=cv, 
                                                    scoring='mean_absolute_error')
rfr_phys_scores = cross_validation.cross_val_score(rfr, phys_features, 
                                                    bandgaps, cv=cv, 
                                                    scoring='mean_absolute_error')

# Let's print out some info about all of these methods performed:
print("Mean absolute error of linear ridge regression on band gap model (Naive):")
print(str(round(abs(mean(scores)), 3)) + " eV")
print("MAE using linear ridge regression and physical features:")
print(str(round(abs(mean(phys_scores)), 3)) + " eV")
print("MAE using random forest regression and naive features:")
print(str(round(abs(mean(rfr_naive_scores)), 3)) + " eV")
print("MAE using random forest regression and physical features:")
print(str(round(abs(mean(rfr_phys_scores)), 3)) + " eV")
print("For comparison, mean absolute error of just guessing avg. band gap is:")
print(str(round(baseline_e, 3)) + " eV")


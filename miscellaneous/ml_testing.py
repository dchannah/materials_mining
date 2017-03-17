import itertools
import numpy as np
from pymatgen import MPRester, periodic_table

"""
Here, we're just messing around with simple machine learning tools from MP.
"""

mp_api = "cVkBCOZrv4TehHfw"

MAX_VECTOR = 100

def get_formula_and_bandgap(chemsys):
    """get_formula_and_bandgap
    Gets the formula and band gap from Materials Project given a chemsys.
    :param chemsys: (String) A binary compound
    :return: [pretty_formula, band_gap]
    """
    stable_results = []
    with MPRester(mp_api) as m:
	ents = m.get_data(chemsys[0] + "-" + chemsys[1], data_type='vasp')
    for e in ents:
	if e['e_above_hull'] < 1e-6:  # Want only thermodynamically stable.
	    stable_results.append([e['pretty_formula'], e['bandgap']])
    return stable_results


def naive_vectorize(comp):
    """naive_vectorize
    Converts a chemical composition into numbers based on atomic fractions.
    :param comp: Composition
    :return: A vector based on the composition's stoichiometry.
    """
    vector = np.zeros(MAX_VECTOR)
    for el in comp:
        frac = comp.get_atomic_fraction(el)
	vector[el.Z - 1] = frac  # Python indexing starts at 0, science doesn't
    return vector


# First, generate a list of all possible binary compounds in the periodic table.
all_binaries = itertools.combinations(periodic_table.all_symbols(), 2)

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

baseline_e = mean(abs(mean(bandgaps)) - bandgaps)
print("Average error", baseline_e)


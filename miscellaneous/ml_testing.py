import itertools
import numpy as np
from pymatgen import MPRester, Composition
from pymatgen.core import periodic_table

"""
Here, we're just messing around with simple machine learning tools from MP.
"""

mp_api = "cVkBCOZrv4TehHfw"
dl = ["pretty_formula", "band_gap", "e_above_hull"]

MAX_VECTOR = 100

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
    print "querying", binary
    queried_data = get_formula_and_bandgap(binary)
    for stable_compound in queried_data:
        materials.append(stable_compound[0])
	bandgaps.append(stable_compound[1])
        naive_features.append(naive_vectorize(stable_compound[0]))

baseline_e = mean(abs(mean(bandgaps)) - bandgaps)
print("Average error", baseline_e)


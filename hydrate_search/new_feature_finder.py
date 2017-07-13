# -*- coding: utf-8 -*-

"""The Feature Finder

A class to do feature finding on compound dictionaries spit out by our
hydrate-finding scripts.

"""

from numpy import mean
from pymatgen import Composition, Element, MPRester

__author__ = "Daniel Hannah"
__email__ = "dchannah@lbl.gov"

class FeatureFinder():
    """
    This class takes a dictionary object from our hydrate-finding codes and
    serves as a wrapper for all of the feature-calculating routines we use to
    do our data mining.

    Attributes:
        comp (dict): The composition of the dry compound in question.
        num_atoms (int): The total number of atoms in the chemical formula.
        mpid (str): The Materials Project ID of the dry compound in question.
        struct (Structure): The atomic structure associated with this MP ID.

    """

    # Lists to categorize cations and anions.  This could be improved.
    working_cations = ["Na", "Li", "H", "K", "Mg", "Ca"]
    poly_cations = ["Si", "Ge", "P", "C"]
    metal_cations = ["Sb", "Zr", "P", "As", "Sn", "Ta", "Y", "Nb", "Ce", "Zn",
                     "Cs", "Bi", "Sc", "Fe", "Dy", "Te", "Sr", "In", "W", "Ga",
                     "Co", "Th", "Ni", "Lu", "Cd", "Cu", "Al", "U", "V", "Yb",
                     "La", "Pr", "Pb", "Gd", "Tm", "Mo", "Rb", "Mn", "Ti", "Pd",
                     "Hg", "Ba", "Hf", "Nd", "Cr", "Be", "Eu", "Br", "B", "Ho",
                     "Er", "Tb", "Sm", "N", "Tl", "Ru", "Ag"]
    all_cations = working_cations + poly_cations + metal_cations
    anions = ["O", "S", "Se", "F", "Cl", "I"]

    def __init__(self, mpid, cmpd_dict):
        """Initializing the feature finder.

        Sets up a FeatureFinder based on a dict from the hydrate finder.

        Args:
            mpid (str): The MP ID of the dry compound in question.
            cmpd_dict (dict): A dictionary output from the hydrate finder.

        """
        comp = Composition(cmpd_dict['formula']).as_dict()
        self.comp = comp
        self.num_atoms = sum([comp[ele] for ele in comp])
        self.mpid = mpid

        with MPRester() as mpr:
            struct = mpr.get_structure_by_material_id(self.mpid)
        self.struct = struct

    def feature_finder(self, feature):
        """A helper method to organize feature derivation calls.

        Args:
            feature (str): The feature we want to calculate.

        Returns:
            The value of the evaluated feature derivation call.
        """
        func_dict = {
            "mean_atomic_radius": self.avg_atomic_radius,
            "sum_atomic_radius": self.sum_atomic_radius,
            "mean_electronegativity": self.avg_x,
            "volume_per_atom": self.vol_per_atom,
            "mean_ionic_radius": self.avg_ionic_radius,
            "sum_ionic_radius": self.ionic_radius_sum,
            "sum_row_number": self.row_number_sum,
            "mean_row_number": self.avg_row_number,
            "sum_group_number": self.group_number_sum,
            "mean_group_number": self.avg_group_number,
            "alkali_metal_ratio": self.alkali_vs_metal,
            "alkali_anion_ratio": self.alkali_vs_anion,
            "metal_anion_ratio": self.metal_vs_anion,
            "coordination_number": self.cation_coordination
        }
        return func_dict[feature]()

    def avg_x(self):
        """Returns the average electronegativity in the compound."""
        sum_x = sum([self.comp[ele] * Element(ele).X for ele in self.comp])
        return sum_x/self.num_atoms

    def sum_atomic_radius(self):
        """Returns the sum of the atomic radii from the formula."""
        sum_ar = sum(self.comp[ele] * Element[ele].atomic_radius for ele in
                     self.comp)
        return sum_ar

    def avg_atomic_radius(self):
        """Returns the average atomic radius from the formula."""
        return self.sum_atomic_radius()/self.num_atoms

    def vol_per_atom(self):
        """Gets a structure from the Materials Project, returns vol/atom."""
        return self.struct.volume/self.struct.num_sites

    def ionic_radius_sum(self):
        """Computes the sum of the average ionic radii in the compound."""
        sum_i_r = sum([self.comp[ele] * Element(ele).average_ionic_radius for
                       ele in self.comp])
        return sum_i_r

    def avg_ionic_radius(self):
        """Returns the average of the average ionic radii in the structure."""
        return self.ionic_radius_sum()/self.num_atoms

    def group_number_sum(self):
        """Returns the sum of the elemental group numbers."""
        sum_group_number = sum([self.comp[ele] * Element(ele).group for ele
                                in self.comp])
        return sum_group_number

    def avg_group_number(self):
        """Returns the average of the elemental group numbers."""
        return self.group_number_sum()/self.num_atoms

    def row_number_sum(self):
        """Returns the sum of the elemental row numbers."""
        r_n_sum = sum([self.comp[ele] * Element(ele).row for ele in self.comp])
        return r_n_sum

    def avg_row_number(self):
        """Returns the average of the elemental row numbers."""
        return self.row_number_sum()/self.num_atoms

    def alkali_vs_metal(self):
        """Ratio of the alkali and metal ionic radii.

        Specifically, this method uses the average of the positive
        oxidation states for the alkali and metal ions.

        Returns:
            The ratio of the alkai and metal ionic radii. 0 if ions not present.
        """
        alkali_radius = 0.0
        metal_radius = 0.0
        for ele in self.comp:
            ox_states = Element(ele).oxidation_states
            if ele in self.working_cations:
                pos_ox_states = [oxs for oxs in ox_states if oxs > 0]
                alkali_radius = mean([Element(ele).ionic_radii[oxs] for oxs
                                      in pos_ox_states if oxs in
                                      Element(ele).ionic_radii])
            if ele in self.metal_cations:
                pos_ox_states = [oxs for oxs in ox_states if oxs > 0]
                metal_radius = mean([Element(ele).ionic_radii[oxs] for oxs
                                     in pos_ox_states if oxs in
                                     Element(ele).ionic_radii])
        if metal_radius == 0.0:
            return 0.0
        else:
            return alkali_radius/metal_radius

    def metal_vs_anion(self):
        """Ratio of the metal radius to the anion radius.

        Specifically, this method uses the average of the positive
        oxidation states for the metal ions and the negative average for anions.

        Returns:
            The ratio of the metal and anion ionic radii.  0 if ions missing.
        """
        metal_radius = 0.0
        anion_radius = 0.0
        for ele in self.comp:
            ox_states = Element(ele).oxidation_states
            if ele in self.metal_cations:
                pos_ox_states = [oxs for oxs in ox_states if oxs > 0]
                metal_radius = mean([Element(ele).ionic_radii[oxs] for oxs
                                     in pos_ox_states if oxs in
                                     Element(ele).ionic_radii])
            if ele in self.anions:
                neg_ox_states = [oxs for oxs in ox_states if oxs < 0]
                anion_radius = mean([Element(ele).ionic_radii[oxs] for oxs in
                                     neg_ox_states if oxs in
                                     Element(ele).ionic_radii])
        if anion_radius == 0.0:
            return 0.0
        else:
            return metal_radius/anion_radius

    def alkali_vs_anion(self):
        """Ratio of the alkali metal radius to the anion radius.

        This method uses the average of the positive oxidation states for the
        alkali and the negative oxidation states for the anion.

        Returns:
            The ratio of the alkali ionic radius and anion ionic radius.
            Returns 0 if ions are missing.
        """
        alkali_radius = 0.0
        anion_radius = 0.0
        for ele in self.comp:
            ox_states = Element(ele).oxidation_states
            if ele in self.working_cations:
                pos_ox_states = [oxs for oxs in ox_states if oxs > 0]
                alkali_radius = mean([Element(ele).ionic_radii[oxs] for oxs
                                      in pos_ox_states if oxs in
                                      Element(ele).ionic_radii])
            if ele in self.anions:
                neg_ox_states = [oxs for oxs in ox_states if oxs < 0]
                anion_radius = mean([Element(ele).ionic_radii[oxs] for oxs in
                                     neg_ox_states if oxs in
                                     Element(ele).ionic_radii])
        if anion_radius == 0.0:
            return 0.0
        else:
            return alkali_radius/anion_radius

    def cation_coordination(self):
        """Gets the average coordination number of cations."""
        b_l = 3.6  # Less than most next-nearest-neighbor distances.
        neighbor_list = []
        for site in self.struct.sites:
            if str(site.specie) in self.all_cations:
                neighbor_list.append(self.struct.get_neighbors(site, b_l))
        avg_coord = mean([len(n_l) for n_l in neighbor_list])
        return avg_coord

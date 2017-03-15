import pd_tools, sys, os, copy, sympy
from fractions import Fraction
from pymatgen import MPRester, Composition, Element, Structure
from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
from pymatgen.phasediagram.maker import PhaseDiagram
from pymatgen.phasediagram.analyzer import PDAnalyzer
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.analysis.structure_matcher import StructureMatcher

"""
This script analyzes the voltages at which cathode materials intercalate working 
ions or convert to a combination of alkali/alkaline earth metals and reduced 
transition metal binaries/ternaries.

Because Materials Project is missing many chalcogenide materials, a complete 
data set necessitates additional local calculations, which we read in here with 
a VaspToComputedEntry drone.  Also, take note that polymporphs are accounted for
in a complete way: Querying for "MnO2", for example, will return many results 
from the MP database.  As a result, intercalation and conversion voltages are 
computed in four different ways:

    * "Maximum": Most stable discharged polymorph + least stable charged
    * "Minimum": Least stable discharged polymorph + most stable charged
    * "Topotactic": Transition metal + anion framework remains unchanged
    * "Ground state": Most stable polymorphs for all states of charge

Note that topotactic matches in the charged and discharged states need not exist
within the stabiliy cutoff.  If no topotactic matches are present in the set of 
structures (the script checks for this), then the topotactic voltages will be 
"None".  The script ultimately writes a CSV file containing the energies and 
IDs of all exact structures used in the voltage calculations.
"""

# Global variables for system stuff
api = "cVkBCOZrv4TehHfw"  # Replace with your MP API.
ehull_cut = {"O": 0.300, "S": 0.365, "Se": 0.450}  # Stability cutoff (eV/atom) 
dl = ["material_id", "spacegroup", "formation_energy_per_atom", "e_above_hull"]

# We also need to define a directory containing supplemental local calculations.
root_dir = "/Users/dchannah/work/tmp/conv_inter_check/calcs_by_sai"

# Chemical systems
monovalents = ["Na", "Li"]
multivalents = ["Mg", "Ca", "Zn"]
working_ions = monovalents + multivalents
tms = ["Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni"]
anions = ["O", "S", "Se"]

def calc_v(starting, resulting):
    """
    Calculates the voltage.
    :param starting: [formula, energy_per_f.u., coefficient] for reactants
    :param resulting: [formula, energy_per_f.u., coefficient] for products
    :return: Voltage for reaction. (float)
    """
    reactant_sum = 0
    product_sum = 0
    n_e = 0
    for pt in resulting:
        f_e_fu = pt[1]
        c = pt[2]
        product = c * f_e_fu
        product_sum += product
    for rt in starting:
        f_e_fu = rt[1]
        c = rt[2]
        product = c * f_e_fu
        reactant_sum += product
        ion_name = str(rt[0])[0:2]
        if ion_name in monovalents:
            n_e = c
        elif ion_name in multivalents:
            n_e = 2 * c
    voltage = -1 * (product_sum - reactant_sum)/n_e
    return voltage


def get_coeffs(dd):
    """get_coeffs
    Returns a dictionary of coefficients based on fractional amounts.
    :param dd: Decomposition dictionary, like {entry: Amount}
    :return: Dictionary of {entry.reduced_formula: coefficient}
    """
    coeff_dict = {}
    for e in dd:
        coeff_dict[e.composition.reduced_composition] = dd[e]
    for r_f in coeff_dict:
        coeff_dict[r_f] = Fraction(str(sympy.nsimplify(coeff_dict[r_f], 
                                   tolerance=1e-4)))
    denominators = [coeff_dict[r_f].denominator for r_f in coeff_dict]
    gcd = max(denominators)
    for r_f in coeff_dict:
        n_a = r_f.reduced_composition.num_atoms
        coeff_dict[r_f] = coeff_dict[r_f] * gcd * (1.0/n_a)
    return coeff_dict, gcd


def get_host_coeff(i_r, c_d, gcd):
    """get_host_coeff
    Gets the coefficient for the "left hand side" of the equation.
    :param i_r: Intercalated entry
    :param c_d: Coeff dict (see get_coeffs)
    :param gcd: Greatest common denominator from get_coeff
    :returns: Coefficient to put in front of WI + host --> etc.
    """
    n_a_host = i_r.composition.reduced_composition.num_atoms
    return gcd/n_a_host


def eles_in_calc(cdir):
    """eles_in_calc
    Checks which elements are in a local calculation.
    :param cdir: Directory containing output files.
    :return: A list of elements.
    """
    s = Structure.from_file(cdir + "/POSCAR")
    return [el for el in s.composition.as_dict()]


def build_complete_pd(f, calc_list):
    """build_complete_pd
    Builds a PD including appropriate local calcs.
    :param f: Formula for which we seek to construct a PD.
    :param calc_list: A list of folders containing local calculations.
    :return: A Phasediagram object.
    """
    entries = pd_tools.get_pruned_chemsys(f)
    csys = [el for el in Composition(f).as_dict()]
    for d in calc_list:
        data_dir = root_dir + "/" + d + "/relax_final"
        if set(eles_in_calc(data_dir)).issubset(set(csys)):
            entries.append(create_entry(data_dir))
    entries_proc = MaterialsProjectCompatibility().process_entries(entries)
    return PhaseDiagram(entries_proc)


def create_entry(cdir):
    """create_entry
    Creates a ComputedStructure entry from a local calculation.
    :param cdir: Directory containing VASP calculation.
    :return: A ComputedStructureEntry object.
    """
    drone = VaspToComputedEntryDrone(inc_structure=True)
    e = drone.assimilate(cdir)
    if e is None:
        print("Missing local calculation for", cdir)
        sys.exit(0)
    else:
        e.entry_id = cdir.rsplit("/", 3)[-2]
        return e


def get_subdir(folder):
    """
    Returns all subdirectories.
    :param folder: Folder for which we want a list of subdirectories.
    :return: List of subdirectories.
    """
    return [name for name in os.listdir(folder) if
            os.path.isdir(os.path.join(folder, name))]


def get_polymorphs(f_str, pda):
    """get_polymorphs
    Gets all the polymorphs with a chemical composition within the cutoff.
    :param f_str: Chemical composition we're interested in.
    :param pda: A PDAnalyzer for the appropriate phase diagram.
    :return: A list of entries within the cutoff.
    """
    comp = Composition(f_str).reduced_composition
    relevant_polymorphs = []
    for e in pda._pd.all_entries:
        if e.composition.reduced_composition == comp:
            ehull = pda.get_decomp_and_e_above_hull(e)[1]
            anion = [el for el in e.composition.as_dict() if el in anions][0]
            if ehull < ehull_cut[anion]:
                form_e = pda._pd.get_form_energy_per_atom(e)
                e.data["formation_energy_per_atom"] = form_e 
                relevant_polymorphs.append(e)
    return relevant_polymorphs


def remove_working_ion(structure):
    """
    Remove all of the working ions from a structure. Looks for any working ion.
    :param structure: A Structure object containing some intercalant.
    :return: The same host structure, but without the working ions.
    """
    s_copy = structure.copy()
    indices_to_remove = []
    for i, s in enumerate(structure.sites):
        if str(s.specie) in working_ions:
            indices_to_remove.append(i)
    s_copy.remove_sites(indices_to_remove)
    return s_copy


def get_voltages(discharged, charged, pda):
    """get_voltages
    Gets relevant voltages for sets of charged/discharged compounds.
    :param discharged: List of intercalated entries.
    :param charged: List of empty host framework entries.
    :param pda: A PDAnalyzer containing the full discharged chemsys.
    :return: V_max, V_min, V_gs, V_topo (optional)
    """
    # We need to sort the entries and get our boundaries.
    disch_e_dict, ch_e_dict = {}, {}
    for e in discharged:
        disch_e_dict[e] = e.data["formation_energy_per_atom"]
    for e in charged:
        ch_e_dict[e] = e.data["formation_energy_per_atom"]
    disch_most_stable = min(disch_e_dict, key=disch_e_dict.get)
    disch_least_stable = max(disch_e_dict, key=disch_e_dict.get)
    ch_most_stable = min(ch_e_dict, key=ch_e_dict.get)
    ch_least_stable = max(ch_e_dict, key=ch_e_dict.get)

    # We always calculate the ground state, max, and min voltages:
    v_max = [get_voltage(disch_most_stable, ch_least_stable, pda),
             {"discharged": disch_most_stable.entry_id,
              "charged": ch_least_stable.entry_id}]
    v_min = [get_voltage(disch_least_stable, ch_most_stable, pda),
            {"discharged": disch_least_stable.entry_id,
             "charged": ch_most_stable.entry_id}]
    v_gs = [get_voltage(disch_most_stable, ch_most_stable, pda),
            {"discharged": disch_most_stable.entry_id,
             "charged": ch_most_stable.entry_id}]

    # We check if there's a topotactic match for disch g.s. for v_topo:
    v_topo = [[None, None, None], None]
    disch_empty = remove_working_ion(disch_most_stable.structure)
    for e in charged:
        sm = StructureMatcher(primitive_cell=False, scale=True,
                              attempt_supercell=True)
        if sm.fit(disch_empty, e.structure):
            v_topo = [get_voltage(disch_most_stable, e, pda),
                      {"discharged": disch_most_stable.entry_id,
                       "charged": e.entry_id}]

    v_dict = {
            "max_intercalation": v_max[0][0],
            "max_conversion": v_max[0][1],
            "max_rxn": v_max[0][2],
            "max_dict": v_max[1],
            "min_intercalation": v_min[0][0],
            "min_conversion": v_min[0][1],
            "min_rxn": v_min[0][2],
            "min_dict": v_min[1],
            "ground_state_intercalation": v_gs[0][0],
            "ground_state_conversion": v_gs[0][1],
            "ground_state_rxn": v_gs[0][2],
            "ground_state_dict": v_gs[1],
            "topotactic_intercalation": v_topo[0][0],
            "topotactic_conversion": v_topo[0][1],
            "topotactic_rxn": v_topo[0][2],
            "topotactic_dict": v_topo[1]
            }

    return v_dict


def get_form_e_per_fu(composition, form_e):
    """
    Calculations (number of atoms in formula) * (formation energy per atom)
    :param composition: A Composition object.
    :param form_e: Formation energy per atom of a compound with composition.
    :return: Formation energy per formula unit.
    """
    sum_atoms = 0
    comp_as_dict = composition.as_dict()
    for key in comp_as_dict:
        sum_atoms += comp_as_dict[key]
    return form_e * sum_atoms


def get_voltage(interc_entry, empty_entry, pda):
    """get_voltage
    Calculates the intercalation and conversion voltage for the entry.
    :param interc_entry: Intercalated entry (e.g. "MgTMO2")
    :param empty_entry: Empty host entry (e.g. "TMO2")
    :param pda: PDAnalyzer containing the full intercalated system.
    :return: V_int, V_conv for the entries in question.
    """
    # First we'll get the decomposition products of the interc entry.
    decomp_dict = get_deco(interc_entry, pda)
    decomp_coeff_dict, rhs_gcd = get_coeffs(decomp_dict)
    interc_entry_coeff = get_host_coeff(interc_entry, decomp_coeff_dict, rhs_gcd)

    # In order to scale the host coefficient correctly, we need N_WI/N_TM 
    num_wi, num_tm = 0, 0
    for ele in interc_entry.composition.as_dict():
        if ele in working_ions:
            num_wi = interc_entry.composition.as_dict()[ele]
        elif ele in tms:
            num_tm = interc_entry.composition.as_dict()[ele]
    host_scaling = num_tm/num_wi
    
    # Get a data triplet for the working ion.
    wi_entry = None
    for wi in working_ions:
        if wi in interc_entry.composition.as_dict():
            wi_entry = pd_tools.get_min_entry(wi)
    wi_f_e = wi_entry.data["formation_energy_per_atom"]
    wi_f_e_fu = get_form_e_per_fu(wi_entry.composition.reduced_composition, 
                                  wi_f_e)
    wi_triplet = [str(wi_entry.composition.reduced_composition), 
                  wi_f_e_fu,
                  interc_entry_coeff]

    # Get a data triplet for the host framework:
    empty_f_e = empty_entry.data["formation_energy_per_atom"]
    empty_f_e_fu = get_form_e_per_fu(empty_entry.composition.reduced_composition,
                                     empty_f_e)
    empty_triplet = [str(empty_entry.composition.reduced_composition),
                     empty_f_e_fu,
                     host_scaling * interc_entry_coeff]

    # Get a data triplet for the intercalated entry:
    interc_f_e = interc_entry.data["formation_energy_per_atom"]
    int_f_e_fu = get_form_e_per_fu(interc_entry.composition.reduced_composition,
                                   interc_f_e)
    interc_triplet = [str(interc_entry.composition.reduced_composition),
                      int_f_e_fu,
                      interc_entry_coeff]

    # Get data triplets for the decomposition products.
    conversion_products = []
    for d in decomp_dict:
        d_f_e = d.data["formation_energy_per_atom"]
        d_f_e_fu = get_form_e_per_fu(d.composition.reduced_composition,
                                     d_f_e)
        d_coeff = decomp_coeff_dict[d.composition.reduced_composition]
        d_triplet = [str(d.composition.reduced_composition),
                     d_f_e_fu,
                     d_coeff]
        conversion_products.append(d_triplet)

    # For now, let's just return the triplets
    int_voltage = calc_v([empty_triplet, wi_triplet], [interc_triplet])
    conv_voltage = calc_v([empty_triplet, wi_triplet], conversion_products)
    return int_voltage, conv_voltage, (interc_triplet, empty_triplet, 
                                       conversion_products)

def get_deco(e, pda):
    """get_deco
    Gets the non-stoichiometric decomposition products of e.
    :param e: Entry we want the decomposition of.
    :param pda: PDAnalzyer containing e.
    """
    e_c = e.composition.reduced_composition
    shifted_entries = [copy.deepcopy(e) for e in pda._pd.all_entries]
    for o_e in shifted_entries:
        if o_e.composition.reduced_composition == e_c:
            o_e.uncorrected_energy += 10.0
    new_pd = PhaseDiagram(shifted_entries)
    new_pda = PDAnalyzer(new_pd)
    decomp = new_pda.get_decomposition(e_c)
    return decomp


def get_deco_old(e, pda):
    """get_deco
    Gets the non-stoichiometric decomposition products of e.
    :param e: Entry we want the decomposition of.
    :param pda: PDAnalzyer containing e.
    """
    # If the entry is stable, we need to lift it off the hull.
    e_c = e.composition.reduced_composition
    if "e_above_hull" in e.data and e.data["e_above_hull"] < 0.008:
        shifted_entries = [copy.deepcopy(e) for e in pda._pd.all_entries]
        for o_e in shifted_entries:
            if o_e.composition.reduced_composition == e_c:
                o_e.uncorrected_energy += 10.0
        new_pd = PhaseDiagram(shifted_entries)
        new_pda = PDAnalyzer(new_pd)
        decomp = new_pda.get_decomposition(e_c)
    elif pda.get_decomp_and_e_above_hull(e)[1] < 0.008:
        shifted_entries = [copy.deepcopy(e) for e in pda._pd.all_entries]
        for o_e in shifted_entries:
            if o_e.composition.reduced_composition == e_c:
                o_e.uncorrected_energy += 10.0
        new_pd = PhaseDiagram(shifted_entries)
        new_pda = PDAnalyzer(new_pd)
        decomp = new_pda.get_decomposition(e_c)
    else:
        decomp = pda.get_decomposition(e_c)
    return decomp


def get_host_coeff(i_r, c_d, gcd):
    """get_host_coeff
    Gets the coefficient for the "left hand side" of the equation.
    :param i_r: Intercalated entry
    :param c_d: Coeff dict (see get_coeffs)
    :param gcd: Greatest common denominator from get_coeff
    :returns: Coefficient to put in front of WI + host --> etc.
    """
    n_a_host = i_r.composition.reduced_composition.num_atoms
    return gcd/n_a_host


def get_coeffs(dd):
    """get_coeffs
    Returns a dictionary of coefficients based on fractional amounts.
    :param dd: Decomposition dictionary, like {entry: Amount}
    :return: Dictionary of {entry.reduced_formula: coefficient}
    """
    coeff_dict = {}
    for e in dd:
        coeff_dict[e.composition.reduced_composition] = dd[e]
    for r_f in coeff_dict:
        coeff_dict[r_f] = Fraction(str(sympy.nsimplify(coeff_dict[r_f], 
                                   tolerance=1e-4)))
    denominators = [coeff_dict[r_f].denominator for r_f in coeff_dict]
    gcd = max(denominators)
    for r_f in coeff_dict:
        n_a = r_f.reduced_composition.num_atoms
        coeff_dict[r_f] = coeff_dict[r_f] * gcd * (1.0/n_a)
    return coeff_dict, gcd

if __name__ == "__main__":
    n_tm = 2
    n_anion = 2 * n_tm
    
    outfile = open('./results.csv', 'w')
    outfile.write("full_compound,working_ion,transition_metal,anion,v_int_max,"  
                  "v_conv_max,max_rxn,max_species,v_int_min,v_conv_min,min_rxn,"
                  "min_species,v_int_gs,v_conv_gs,gs_rxn,gs_species,v_int_topo,"
                  "v_conv_topo,v_conv_rxn,v_conv_species\n")
    
    calclist = get_subdir(root_dir)
    for ion in working_ions:
        for tm in tms:
            for anion in anions:
                interc_formula = ion + tm + str(n_tm) + anion + str(n_anion)
                print("About to run", interc_formula)
                empty_formula = tm + str(n_tm) + anion + str(n_anion)
                pd_full = build_complete_pd(interc_formula, calclist)
                pda_full = PDAnalyzer(pd_full) 
                structs_interc = get_polymorphs(interc_formula, pda_full)
                structs_empty = get_polymorphs(empty_formula, pda_full)
                voltage_dict = get_voltages(structs_interc, structs_empty, pda_full)
                outfile.write(str(interc_formula) + "," + str(ion) + "," + str(tm)\
                        + "," + str(anion) + "," + 
                        str(voltage_dict["max_intercalation"]) + "," + 
                        str(voltage_dict["max_conversion"]) + "," +
                        str(voltage_dict["max_rxn"]) + "," + 
                        str(voltage_dict["max_dict"]) + "," + 
                        str(voltage_dict["min_intercalation"]) + "," + 
                        str(voltage_dict["min_conversion"]) + "," + 
                        str(voltage_dict["min_rxn"]) + "," + 
                        str(voltage_dict["min_dict"]) + "," + 
                        str(voltage_dict["ground_state_intercalation"]) + "," + 
                        str(voltage_dict["ground_state_conversion"]) + "," +
                        str(voltage_dict["ground_state_rxn"]) + "," + 
                        str(voltage_dict["ground_state_dict"]) + "," + 
                        str(voltage_dict["topotactic_intercalation"]) + "," + 
                        str(voltage_dict["topotactic_conversion"]) + "," + 
                        str(voltage_dict["topotactic_rxn"]) + "," + 
                        str(voltage_dict["topotactic_dict"]) + "\n")

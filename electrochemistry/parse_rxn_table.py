from pymatgen import Composition
from formula_to_tex import latex_from_formula, working_ions, tms


def build_coeff_dict(l):
    """
    Builds the appropriate dictionary from a list of decomp products.
    :param l: Doubles put out from my v_calc script.
    :return: {formula: decomp coefficient} for everything in l
    """
    coeff_dict = {}
    for triplet in l:
        coeff_dict[triplet[0]] = triplet[2]
    return coeff_dict


def parse_line(l):
    """
    Parses a line from the output CSV file from v_calc.
    :param l: ([interc, e, amount], [ch, e, amt], [[conv_1, e, amt]...[conv_N...
    :return: A chemical reaction based on this line from input file (LaTeX)
    """
    # Convert the line to Python list objects.
    list_of_lists = eval(l)

    # Find the working ion
    interc_formula_dict = Composition(list_of_lists[0][0]).as_dict()
    working_ion = [el for el in interc_formula_dict if el in working_ions][0]
    working_ion_coeff = list_of_lists[0][2]

    # Read the host formula from the list objects.
    host_formula = list_of_lists[1][0]
    host_coeff = list_of_lists[1][2]

    # Build the intercalation in LaTeX
    tm = [el for el in interc_formula_dict if el in tms][0]
    if int(interc_formula_dict[tm]) == 2:
        int_rxn_str = working_ion + " + $\\frac{1}{2}$" + \
                      latex_from_formula(host_formula)
    else:
        int_rxn_str = working_ion + " + " + latex_from_formula(host_formula)
    int_rxn_str += " $\\rightarrow$ " + latex_from_formula(list_of_lists[0][0])

    # Get the list of decomp. product triplets and convert to a dictionary.
    list_of_decomp_products = list_of_lists[2]
    decomp_dict = build_coeff_dict(list_of_decomp_products)

    # Build the conversion reaction in LaTeX
    if working_ion_coeff != 1:
        if working_ion_coeff == 0.5:
            conv_rxn_str = "$\\frac{1}{2}$" + working_ion + " + "
        else:
            conv_rxn_str = str(int(working_ion_coeff)) + working_ion + " + "
    else:
        conv_rxn_str = working_ion + " + "
    if host_coeff == 0.5:
        conv_rxn_str += "$\\frac{1}{2}$" + latex_from_formula(host_formula)
    else:
        conv_rxn_str += str(int(host_coeff)) + latex_from_formula(host_formula)
    conv_rxn_str += " $\\rightarrow$ "
    first_product = True
    for conv_product in decomp_dict:
        if not first_product:
            conv_rxn_str += " + "
        tex_formula = latex_from_formula(conv_product)
        if decomp_dict[conv_product] == 0.5:
            conv_rxn_str += "$\\frac{1}{2}$" + tex_formula
        elif decomp_dict[conv_product] != 1:
            conv_rxn_str += str(int(decomp_dict[conv_product])) + tex_formula
        else:
            conv_rxn_str += tex_formula
        if first_product:
            first_product = False

    tex_str = int_rxn_str + " & " + conv_rxn_str + " \\\\"
    return tex_str


def read_results_dict(filename):
    with open(filename, 'r') as f:
        for l in f.readlines():
            print(parse_line(l))

if __name__ == "__main__":
    test_file = "./one_electron_reactions.txt"
    read_results_dict(test_file)


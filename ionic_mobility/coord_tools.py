import sys
import numpy as np
from pymatgen import Structure


def get_energy(cn, x, planar=False):
    """get_energy
    Gets the energy associated with a particular BL and CN. 
    :param cn: The coordination number of the site. 
    :param x: The bond distance at that site.
    :param planar: Describes whether coordination is planar.
    :return: Energy interpolated from the curve.
    """
    if cn == "4" and planar:
        cn = "4sq"
    elif cn == "4" and not planar:
        cn = "4tet"
    else:
        cn = str(cn)
    cn_energy_dict = {
            "2": (-68.6293884644*(x**5) + 656.9448984491*(x**4) - 2547.6671761795*(x**3) + 5012.2963554218*(x**2) - 5004.1902727195*x + 2015.9081144944),
            "3": (-34.9312858121*(x**5) + 376.7504178486*(x**4) - 1639.9994652078*(x**3) + 3606.7810071703*(x**2) - 4005.8593944802*x + 1783.0421752158),
            "4tet": (-25.3423639852*(x**5) + 288.7738787470*(x**4) - 1327.4907307968*(x**3) + 3081.2192643508*(x**2) - 3608.0587468723*x + 1689.7838649931),
            "4sq": (-6.3642389833*(x**5) + 83.0770384949*(x**4) - 436.8216229413*(x**3) + 1155.6813299862*(x**2) - 1530.2182605825*x + 795.0068003544),
            "5": (-19.7300654532*(x**5) + 235.9012586699*(x**4) - 1139.2117393492*(x**3) + 2781.5508184700*(x**2) - 3431.7317287412*x + 1696.4278974512),
            "6": (-16.3481342917*(x**5) + 200.5759354262*(x**4) - 993.4589156796*(x**3) + 2485.6204381302*(x**2) - 3137.2859206340*x + 1581.6198097747),
            "8": (-6.6191602186*(x**5) + 92.6879849303*(x**4) - 522.2274246236*(x**3) + 1480.7902801777*(x**2) - 2108.2467380805*x + 1190.3906582120),
            }
    return cn_energy_dict[cn]


def get_average_bl(struct, idx, cn):
    """get_average_bl
    Gets the avergae bond distance for a coordination number at index.
    :param struct: Structure
    :param idx: Index we want the bond distance from.
    :param cn: Coordination number of the atom in question.
    :return: Average bond distance of the top CN neighbors.
    """
    cn = int(cn)
    distances = []
    for i, s in enumerate(struct.sites):
        if i != idx and str(s.specie) in ["O2-", "S2-", "Se2-"]:
            distances.append(struct.get_distance(i, idx))
    distances = np.sort(np.array(distances))
    coordinated = distances[:cn]
    return np.average(coordinated)

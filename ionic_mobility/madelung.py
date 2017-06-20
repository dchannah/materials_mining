import sys
import math
import numpy as np
from pymatgen import Structure, Site, Specie
from coord_tools import get_average_bl, get_energy
from oxi_filter import compute_tm_ox, tms


def create_vac(supercell, orig_vac, orig_lat):
    """create_vac
    Creates an equivalent vacancy in a supercell.
    :param supercell: A Pymatgen Structure representation of the supercell.
    :param orig_vac: Original vacancy coordinates.
    :param orig_lat: Original lattice - Pymatgen lattice object.
    :return: None
    """
    new_vacancy_index = scale_coord_index(supercell, orig_vac, orig_lat)
    if new_vacancy_index is not None:
        supercell.remove_sites([new_vacancy_index])
        # supercell.replace(new_vacancy_index, Specie.from_string("Ca2+"))
    else:
        print "Nothing got removed!"
        sys.exit(0)
    return


def shift_to_activated(supercell, mg_coords, act_coords, orig_lat):
    """shift_to_activated
    Shifts the Mg atom to the activated state coordinates in a supercell.
    :param supercell: A pymatgen structure representation of the supercell.
    :param mg_coords: Mg coordinates in the original cell.
    :param act_coords: Activated state coordinates in the original cell.
    :param orig_lat: Original lattice - Pymatgen lattice object
    :return: None
    """
    mg_idx = scale_coord_index(supercell, mg_coords, orig_lat)
    act_coords = scale_coords(supercell, act_coords, orig_lat)
    supercell.replace(mg_idx, supercell.sites[mg_idx].specie, coords=act_coords)
    return


def scale_coords(supercell, coords, orig_lat):
    """scale_coords.
    Scales coordinates in a small cell to a bigger one.
    :param supercell: A Pymatgen Structure representation of the supercell.
    :param coords: Coordinates we seek to scale.
    :param orig_lat: Original lattice - Pymatgen lattice object.
    :return: [a_new, b_new, c_new]
    """

    # Define a method to cleanly compute RMS.
    def rms(x, axis=None):
        return np.sqrt(np.mean(x**2, axis=axis))

    a_new = coords[0] * rms(orig_lat.matrix[0])/rms(supercell.lattice.matrix[0])
    b_new = coords[1] * rms(orig_lat.matrix[1])/rms(supercell.lattice.matrix[1])
    c_new = coords[2] * rms(orig_lat.matrix[2])/rms(supercell.lattice.matrix[2])
    return [a_new, b_new, c_new]


def scale_coord_index(supercell, coords, orig_lat):
    """scale_coord_index
    Finds the index of the equivalent (scaled) site in the supercell.
    :param supercell: A Pymatgen Structure representation of the supercell.
    :param coords: Coordinates we seek to scale.
    :param orig_lat: Original lattice - Pymatgen lattice object.
    :return: A Site index
    """
    tol = 1e-4
    new_coords = scale_coords(supercell, coords, orig_lat) 
    idx = None
    for i, s in enumerate(supercell.sites):
        diffx = abs(np.array(s._fcoords)[0] - np.array(new_coords)[0])
        diffy = abs(np.array(s._fcoords)[1] - np.array(new_coords)[1])
        diffz = abs(np.array(s._fcoords)[2] - np.array(new_coords)[2])
        if diffx < tol and diffy < tol and diffz < tol:
            idx = i
            tol = min(diffx, diffy, diffz)
    return idx


def madelung_energy(struct, idx, ref_length):
    """madelung_energy
    Computes the Madelung energy of a specified site.
    :param struct: Pymatgen Structure object.
    :param idx: A site index for which we seek the Madelung energy.
    :param ref_length: A reference length parameter (usually 1)
    :return: Madelung site energy in eV.
    """
    m = 0
    for i, s in enumerate(struct.sites):
        if i != idx:
            z_i = float(s.specie.oxi_state)
            r_ij = struct.get_distance(idx, i)
            m += z_i / (r_ij/ref_length)
    charge = float(struct.sites[idx].specie.oxi_state)
    electron_charge = 1.6022e-19
    four_pi_e0 = 1.112e-10
    ang_to_meter = 1e-10
    r_l_meters = ref_length * ang_to_meter
    solid_eps = 10.0  # Based on https://srd.nist.gov/JPCRD/jpcrd28.pdf
    prefactor_joules = (charge * electron_charge**2)/(four_pi_e0 * solid_eps * r_l_meters)
    e_joules = prefactor_joules * m
    joule_to_ev = 6.242e+18
    e = e_joules * joule_to_ev
    return e


def apply_oxi_states(struct, auto=True):
    """apply_oxi_states
    Applies oxidations states assuming Mg is +2 and anion is -2.
    :param struct: A Structure containing (Mg, TM, anion)
    :return: None
    """
    if auto:
        tm = [el for el in struct.composition.as_dict() if el in tms][0]
        tm_ox = compute_tm_ox(struct.composition.as_dict())
        print "Assigning TM ox as ", tm_ox
        ox_dict = {"Mg": 2, "O": -2, "P": 5, "Si": 4, tm: tm_ox}
        struct.add_oxidation_state_by_element(ox_dict)
    else:
        #user_string = {"Mg": 2, "O": -2, "Ti": 3.5, "P": 5}
        user_string = {"Mg": 2, "O": -2}
        struct.add_oxidation_state_by_element(user_string)
    return



def read_coord_file(filename, dlimit_str="  "):
    """read_coord_file
    Reads in coordinate file from filename.
    :param filename: File containing (mg coords, vac coords, [act coords...])
    :param dlimit_str: Delimiter for reading in text.
    """
    all_coords = np.loadtxt(fname=filename, delimiter=dlimit_str)
    ion_coords = all_coords[0]
    vacancy_coords = all_coords[1]
    act_coords_x = []
    act_coords_y = []
    act_coords_z = []
    for i in range(2,len(all_coords)):
        act_coords_x.append(all_coords[i][0])
        act_coords_y.append(all_coords[i][1])
        act_coords_z.append(all_coords[i][2])
    activated_coords = [np.average(act_coords_x), np.average(act_coords_y), 
                        np.average(act_coords_z)]
    return ion_coords, vacancy_coords, activated_coords


if __name__ == "__main__":
    # Read in the structure from user input
    s = Structure.from_file(sys.argv[1])
    
    # Also read in coordination for stable and activated.
    stable_cn = sys.argv[2]
    activated_cn = sys.argv[3]

    # Need to check if user has specified a planar activated state.
    planar_activated = False
    if len(sys.argv) > 4:
        planar_activated = sys.argv[4]

    # Read in coordinates to scale up diffuser and vacancy indices.
    fname = "coord_file.txt"
    mg_coords, vac_coords, act_coords = read_coord_file(fname)
    apply_oxi_states(s, True)
    #print "Act coords are", act_coords

    # Get the converged Madelung site energy difference.
    deltaE_prev = 99999
    for n in range(1, 100):
        # Make copy and supercell
        scell = s.copy()
        scell.make_supercell([n, n, n])

        # Create vacancy
        create_vac(scell, vac_coords, s.lattice)
        # scell.to(fmt='poscar', filename="with_vac.vasp") 
        # Get Mg index and compute stable site energy.
        mg_idx = scale_coord_index(scell, mg_coords, s.lattice)
        e_stable = madelung_energy(scell, mg_idx, 1)
        
        # On our first run, with the vacancy created, we get stable site cluster E.
        if n == 1:
            scell.to(fmt='poscar', filename="stable.vasp")
            bond_length = get_average_bl(scell, mg_idx, stable_cn)
            print "Stable bond length is", bond_length
            stable_cluster_e = get_energy(stable_cn, bond_length)
        
        # Move the Mg at mg_idx to activated state and re-compute energy.
        shift_to_activated(scell, mg_coords, act_coords, s.lattice)
        # scell.to(fmt='poscar', filename=(str(n) + ".vasp"))
        e_act = madelung_energy(scell, mg_idx, 1)
        
        # Now that Mg is in the activated site, get the cluster E for this site.
        if n == 1:
            act_bond_length = get_average_bl(scell, mg_idx, activated_cn)
            print "Activated site bond length is", act_bond_length
            act_cluster_e = get_energy(activated_cn, act_bond_length, planar_activated)
        
        # Print when converged and break
        e_tol = 0.01
        deltaE = e_act - e_stable
        print n, e_stable, e_act, (e_act - e_stable)
        if np.abs(deltaE - deltaE_prev) < e_tol:
            # print n, e_stable, e_act, (e_act - e_stable)
            break
        # scell.to(fmt='poscar', filename=(str(n) + "_act.vasp"))
        deltaE_prev = deltaE

    # Now compute the final barrier height.
    print "final deltaE", deltaE
    print "Activated cluster energy is ", act_cluster_e
    print "Stable cluster energy is ", stable_cluster_e
    delta_cluster = act_cluster_e - stable_cluster_e
    print "delta cluster is", delta_cluster
    deltaE += delta_cluster
    print("Barrier height (eV): " + str(deltaE))


import numpy as np
from scipy.spatial import ConvexHull
from pymatgen import Molecule


def read_xyz(xyz_file):
    with open(xyz_file, 'r') as f:
        coords = np.genfromtxt(f, usecols=range(1,4), skip_header=2)
    with open(xyz_file, 'r') as f2:
        species = np.loadtxt(f2, skiprows=2, usecols=[0], dtype=str)
    return coords, species


def sph2cart(r, el,az):
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    return x, y, z


def cart2sph(x, y, z):
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    el = np.arctan2(z, hxy)
    az = np.arctan2(y, x)
    return r, el, az


def cart2sph_all(pts):
    return np.array([cart2sph(x,y,z) for x,y,z in pts])


def sph2cart_all(relevaz):
    return np.array([sph2cart(r, elev, az) for r, elev, az, in relevaz])


def get_surface_indices_1(points):
    hull = ConvexHull(points)
    return hull.vertices


def get_surface_indices_2(xyz_file, bond_distance):
    mol = Molecule.from_file(xyz_file)
    bulk_coordination = 4
    surface_i = []
    for i in range(len(mol.sites)):
        neighbors = mol.get_neighbors(mol.sites[i], bond_distance)
        if len(neighbors) < bulk_coordination:
            surface_i.append(i)
    return surface_i


def get_ligand_coords(cart_coords, surface_i, elements):
    sph_coords = cart2sph_all(cart_coords)
    ligand_sph_coords = []
    ligand_identities = []
    for i in range(len(sph_coords)):
        if i in surface_i:
            ligand_name = "H"
            ligand_radius = 0.53
            surface_species = ligand_name + elements[i]
            ligand_identities.append(surface_species)
            ligand_bond_length = 2.5
            r_ligand = sph_coords[i][0] + ligand_bond_length
            elev_ligand = sph_coords[i][1]
            az_ligand = sph_coords[i][2]
            ligand_sph_coords.append((r_ligand, elev_ligand, az_ligand))
    ligand_cart_coords = sph2cart_all(np.array(ligand_sph_coords))
    return ligand_cart_coords, ligand_identities


def write_passivated_particle(crystal_coords, ligand_coords, elements,
                              lig_species):
    output = open("./surface.xyz", 'w')
    output.write(str(len(crystal_coords) + len(ligand_coords)) + "\n\n")
    for i in range(len(crystal_coords)):
        x, y, z = crystal_coords[i][0], crystal_coords[i][1], \
                  crystal_coords[i][2]
        output.write(elements[i] + " " + str(x) + " " + str(y) + " " + str(z) +
                     "\n")
    for i in range(len(ligand_coords)):
        x, y, z = ligand_coords[i][0], ligand_coords[i][1], ligand_coords[i][2]
        output.write(lig_species[i] + " " + str(x) + " " + str(y) + " " +
                     str(z) + "\n")

xyz_file = "./cdse.xyz"
coordinates, elements = read_xyz(xyz_file)
surface_indices = get_surface_indices_2(xyz_file, 2.8)
# surface_indices = get_surface_indices_1(coordinates)
ligand_coords, ligand_elements = get_ligand_coords(coordinates,
                                                   surface_indices, elements)
write_passivated_particle(coordinates, ligand_coords, elements,
                          ligand_elements)


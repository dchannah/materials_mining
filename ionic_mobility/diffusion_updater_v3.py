import os
from pymatgen import Structure, Site
from pymatgen.io.vasp.outputs import Oszicar, Poscar
from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar
from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
# import pymatgen.io.smartio as writer
import math
import numpy as np

__author__ = 'dchannah,saigautam,pcanepa'
__version__ = "0.2"
__maintainer__ = "dchannah"
__email__ = "dchannah@lbl.gov, gautam91@mit.edu, pcanepa@lbl.gov"
__date__ = "Oct 20, 2015"


class UpdaterNEB:
    def __init__(self, mp_id, working_ion, pair_index):
        self._mp_id = mp_id
        self._working_ion = working_ion
        if len(pair_index) == 1:
            self._pair_index = "00" + pair_index
        elif len(pair_index) == 2:
            self._pair_index = "0" + pair_index
        else:
            self._pair_index = pair_index

    def checkint(self, s):
        """
        @type s: Object to check for integer.

        @return: Returns a Boolean describing whether s is integer.
        """
        try:
            int(s)
            return True
        except ValueError:
            return False

    def checkstart(self, dirslist, root_path='./'):
        """
        Checks if there is a directory numbered 00 and has a POSCAR inside
        """

        if '00' in dirslist and os.path.isfile(root_path + '00/POSCAR') and os.path.isfile(root_path + '00/OSZICAR'):
            return True
        else:
            return False

    def checkend(self, dirslist, root_path='./'):
        """
        Known bug: Both CONTCARS and POSCARS are present in images. Check may not be rigorous.
        """
        if os.path.isfile(root_path+dirslist[-1] + '/POSCAR') and os.path.isfile(root_path+dirslist[-1] + '/OSZICAR'):
            return True
        else:
            return False

    def listdirs(self, folder):
        return [d for d in os.listdir(folder) if os.path.isdir(os.path.join(folder,d))]

    def dirlist(self, root_path='./'):
        """
        Returns a list of directories containing all the images and the start and end points of NEB calculations.

        root_path = Path to the directory containing the results of the NEB run.

        Please make sure you have the naming sequence for the images as XY where X and Y are integers.
        Let 00 represent the starting point and the last integer (09 for example in a system of 8 images)
        represent the end point.
        """
        dirs = self.listdirs(root_path)
        dirs_list = []
        for i in range(0,len(dirs)):
            if self.checkint(dirs[i]) and len(dirs[i])==2:
                dirs_list.append(dirs[i])
        return dirs_list

    def sort_dirs(self, dirslist):
        dirs_new = []
        for i in dirslist:
            dirs_new.append(int(i))
            dirs_new.sort()
        for j in range(0,len(dirslist)):
            if len(str(dirs_new[j])) == 1:
                dirslist[j] = '0'+str(dirs_new[j])
            elif len(str(dirs_new[j])) == 2:
                dirslist[j] = str(dirs_new[j])
            else:
                assert 0 == 1, "Check directory numbering or naming. It should be of the form XY where XY is integer"
        return dirslist

    def get_diffuser_index(self, structures, cation):
        diff_ion_index = 0
        max_path_length = 0
        for i in range(0, len(structures[0].sites)):
            if str(structures[0].sites[i].specie) == cation:
                path_length = self.get_path_length(structures, i)
                if path_length > max_path_length:
                    max_path_length = path_length
                    diff_ion_index = i
        return diff_ion_index

    def get_diff_sp(self, poscars, atommap, user_sp='', user_pos=None):
        """
        Returns diffusing species and it's position

        Known bug: Assumes the diffusing species is dilute, i.e., it's the kind of element that's present in least
        quantity in the POSCAR. This may not be a rigorous check. Please pass in user_sp if results are contradicting.
        Might be difficult to track the diffusing species in dilute vacancy limit.

        user_sp: If you already know the diffusing species or want to fix a particular element, pass it here.

        user_pos: If there are multiple diffusing ions in the lattice, use this variable to choose one for which
                the bond distances should be calculated. Pass an integer, string or a single level list. One can
                extract user_pos from the dictionary of user_pos and atommap
        """
        poscar = poscars[0]
        no_atoms = poscar.natoms
        dict_atoms = {}
        element = poscar.site_symbols
        for i in range(0,len(element)):
            dict_atoms.update({element[i]:int(no_atoms[i])})
            no_atoms[i] = int(no_atoms[i])
        no_atoms.sort()
        if user_sp != '':
            sp = user_sp
        elif no_atoms[0] != no_atoms[1] and no_atoms[0] == 1:
            sp = dict_atoms.keys()[dict_atoms.values().index(no_atoms[0])]
        else:
            assert 0 == 1, "There is at least 1 more element in the unit cell with the same number of atoms as the " \
                           "diffusing species. Specify user_sp and user_pos if this is the case."

        pos = []
        if user_pos is not None:
            if type(user_pos) is list:
                posn = user_pos
                for d in posn:
                    pos.append(int(d))
            elif self.checkint(user_pos):
                pos.append(int(user_pos))
            else:
                assert 0==1, "user_pos must be an int, str or list"
        else:
            for j in atommap:
                s = ''
                posn = s
                for k in range(0,len(j)):
                    if not self.checkint(j[k]):
                        s += j[k]
                    else:
                        posn += j[k]
                if sp == s:
                    pos.append(int(posn))
        if sp+str(pos[0]) not in atommap:
            assert 0 == 1, "Values supplied in user_sp and user_pos dont match Poscar values. Check input"
        else:
            diff_sp = {sp: pos}
        return diff_sp

    def get_actual_diff_sp(self, poscars, all_cart_list, atommap, user_sp=' ', user_pos=None):
        """
        Caveat: Returns the atom that moves the longest. It shouldn't matter in coupled diffusion. But if you do want
        two different path lengths, pass the user_sp and user_pos seperately

        In case of doubt, always pass the user_sp and user_pos explicitly
        """
        if user_sp == ' ' and user_pos is None:
            prev_dist = 0.0
            prev_posn = 0
            prev_s = ' '
            for j in atommap:
                s = ''
                posn = s
                for k in range(0,len(j)):
                    if not self.checkint(j[k]):
                        s += j[k]
                    else:
                        posn += j[k]
                posn = int(posn)
                cart = all_cart_list[0][posn-1]
                cart1 = all_cart_list[-1][posn-1]
                dist = math.sqrt((cart[0]-cart1[0])**2+(cart[1]-cart1[1])**2+(cart[2]-cart1[2])**2)
                if dist > prev_dist:
                    prev_dist = dist
                    prev_posn = posn
                    prev_s = s
        else:
            prev_s = user_sp
            prev_posn = user_pos
        diff_sp = self.get_diff_sp(poscars, atommap, prev_s, prev_posn)
        return diff_sp

    def get_cart_images(self, dirslist):
        """
        Gets cartesian coordinates of all images in the NEB as a list.

        dirlist= List of directories that contain the POSCARs or the CONTCARs

        """
        cart_all = []
        poscar = self.get_all_poscars(dirslist)
        for i in range(0,len(dirslist)):
            cart_all.append(poscar[i].structure.lattice.get_cartesian_coords(poscar[i].structure.frac_coords))
        return cart_all

    def get_frac_images(self, dirslist):
        """
        Gets fractional coordinates of all images in the NEB as a list.

        dirlist= List of directories that contain the POSCARs or the CONTCARs

        """
        frac_all = []
        poscar = self.get_all_poscars(dirslist)
        for i in range(0,len(dirslist)):
            frac_all.append(poscar[i].structure.frac_coords)
        return frac_all

    def get_atommap(self, dirslist, poscars_list):
        atommap = []
        if self.checkstart(dirslist):
            poscar = poscars_list[0]
            no_atoms = poscar.natoms
            element = poscar.site_symbols
            count = 1
            for k in range(0,len(element)):
                if k != 0:
                    count += int(no_atoms[k-1])
                for l in range(0,int(no_atoms[k])):
                    atommap.append(element[k]+str(l+count))
        else:
            assert 0 == 1, "There is no directory number 00. Check your directory lists"
        return atommap

    def get_cart_from_frac(self, coords, lattice):
        """
        Gets cartesian coordinates for a set of fractional coordinates and lattice matrix passed.

        Very trivial function. Ensure that the order of fractional coordinates and the lattice match. No checks made.

        """
        cart_coords = []
        for i in range(0,len(coords)):
            cart_coords.append(coords[i]*(lattice[0][i]+lattice[1][i]+lattice[2][i]))
        return cart_coords

    def get_lattice_constant(self, poscars):
        """
        Function to get lattice constants. Uses pymatgen.core.lattice; returns a list of [x,y,z]
        """
        if type(poscars) is list:
            poscar = poscars[0]
        else:
            poscar = poscars
        return poscar.structure.lattice._lengths

    def get_lattice(self, poscars):
        """
        Function to get lattice constant matrix. Returns a 2d list if one POSCAR is passed. If a list of
        POSCARs are passed, then returns a 3D-list
        """
        if type(poscars) is not list:
            return poscars.structure.lattice._matrix
        else:
            matrix = []
            for i in poscars:
                matrix.append(i.structure.lattice._matrix)
            return matrix

    def get_cn(self, poscars, diff_sp):
        """
        A function to compute the coordinations numbers not rounded

        Output: cn
        """
        n = int(diff_sp.values()[0][0])-1
        cn=[]
        for i in range(0,len(poscars)):
            coord_no = VoronoiCoordFinder(poscars[i].structure)
            cn.append(coord_no.get_coordination_number(n))
        return cn

    def vasp_neb_to_pymat(self, dirs):
        """
        Takes structures from NEB images and returns an array of pymatgen Structure objects.
        Directory list must be in numerical order along NEB path (i.e. 00, 01, ..., N)
        """
        s_array = []
        for i in range(0, len(dirs)):
            if i == 0 or i == (len(dirs) - 1):
                s_array.append(Structure.from_file(dirs[i] + "/POSCAR"))
            else:
                s_array.append(Structure.from_file(dirs[i] + "/CONTCAR"))
        return s_array

    def get_all_poscars(self, dirslist,root_path='./'):
        poscar = []
        if self.checkstart(dirslist) and self.checkend(dirslist):
            for i in range(0,len(dirslist)):
                if i != 0 and i != len(dirslist)-1:
                    path = os.path.join(root_path, dirslist[i], "CONTCAR")
                else:
                    path = os.path.join(root_path, dirslist[i], "POSCAR")
                poscar.append(Poscar.from_file(path))
        else:
            assert 0 == 1, "Check your start and end directories. POSCAR or OSZICAR may be missing"
        return poscar

    def get_empty_lattice(self, structures, cation):
        starting_struct = structures[0].copy()
        indices_to_remove = []
        for s in range(0, len(starting_struct.sites)):
            if str(starting_struct.sites[s].specie) == cation:
                indices_to_remove.append(s)
        starting_struct.remove_sites(indices_to_remove)
        return starting_struct

    def get_start_end_sites(self, structures, index):
        starting_struct, ending_struct = structures[0], structures[len(structures)-1]
        return starting_struct.sites[index], ending_struct.sites[index]

    def get_cation_concentration(self, structures, cation):
        pymat_structure = structures[0]
        num_diffusers = float(len(pymat_structure.indices_from_symbol(cation)))
        tm_ions = ["Mn", "V", "Fe", "Ti", "Cr", "Co", "Mo", "Ni", "Cu", "Bi", "Sb", "Ta", "W", "Re", "Sn"]
        num_tm = 0
        for tm in tm_ions:
            num_of_this_one = float(len(pymat_structure.indices_from_symbol(tm)))
            num_tm += num_of_this_one
        return num_diffusers/num_tm

    def get_path_sites(self, structures, index):
        path_sites = []
        for i in range(0, len(structures)):
            path_sites.append(structures[i].sites[index].as_dict())
        return path_sites

    def get_neighbors(self, poscars, atommap, voronoi_sites, tol = 1e-6):
        """
        Gets the indices of the neighboring atoms near the diffusing species (or any other species whose voronoi_sites
        are known.

        The indices are returned as a list for each poscar passed.

        The indices have the same numerical value as atommap/user_pos
        """
        neighbors = []
        for l in range(0,len(poscars)):
            array = []
            count = 'False'
            for k in range(0,len(voronoi_sites[l])):
                for j in range(0,len(atommap)):
                    for i in range(0,3):
                        if abs(voronoi_sites[l][k]._fcoords[i]-poscars[l].structure.frac_coords[j][i]) < tol:
                            count = 'True'
                        else:
                            if abs(voronoi_sites[l][k]._fcoords[i]-(poscars[l].structure.frac_coords[j][i]-1.0)) < tol:
                                poscars[l].structure.frac_coords[j][i] -= 1.0
                                count = 'True'
                            elif abs(voronoi_sites[l][k]._fcoords[i]-(poscars[l].structure.frac_coords[j][i]+1.0)) \
                                    < tol:
                                poscars[l].structure.frac_coords[j][i] += 1.0
                                count = 'True'
                            else:
                                count = 'False'
                                break
                    if count == 'True':
                        array.append(int(j)+1)
            neighbors.append(array)
        return neighbors

    def get_voronoi(self, poscars, diff_sp):
        """
        Known bug, solid_angle_tol = 0.5 works for Mn2O4 spinels. Needs to be checked for others.
        solid_angle_tol is the tolerance on the solid angle, bigger the angle- closer the atoms.
        """
        solid_angle_tol = 0.5
        sites = []
        n = int(diff_sp.values()[0][0])-1
        for i in range(0, len(poscars)):
            voronoi = VoronoiCoordFinder(poscars[i].structure)
            sites.append(voronoi.get_coordinated_sites(n, solid_angle_tol))
        return sites

    def get_nn_bond_length(self, poscars, neighbors, diff_sp):
        """
        Returns a list of bond lengths required for bond valence method on each image.
        Alert : 2D list is returned
        """
        bond_lengths = []
        n = int(diff_sp.values()[0][0])-1
        for j in range(0, len(poscars)):
            dist = []
            for i in range(0, len(neighbors[j])):
                sp = int(neighbors[j][i])-1
                dist.append(poscars[j].structure.get_distance(n, sp))
            bond_lengths.append(dist)
        return bond_lengths

    def get_all_oszicars(self, dirslist, root_path='./'):
        oszicar = []
        if self.checkstart(dirslist) and self.checkend(dirslist):
            for i in range(0, len(dirslist)):
                path = os.path.join(root_path, dirslist[i], "OSZICAR")
                oszicar.append(Oszicar(path))
        else:
            assert 0 == 1, "Check your start and end directories. OSZICAR may be missing"
        return oszicar

    def get_energies_images(self, dirslist):
        energies = []
        oszicar = self.get_all_oszicars(dirslist)
        for i in range(0, len(dirslist)):
            energies.append(oszicar[i].ionic_steps[-1]["E0"])
        return energies

    def get_scaled_energies(self, energies_images):
        energy_endpoint_mean = (energies_images[0] + energies_images[-1])*0.5
        scaled_energies = np.asarray(energies_images)
        scaled_energies = (scaled_energies-energy_endpoint_mean) * 1000
        return scaled_energies

    def max_energy(self, energy_array):
        max_energy = float(energy_array[0])
        for i in range(1, len(energy_array)):
            if float(energy_array[i]) > max_energy:
                max_energy = energy_array[i]
        return max_energy

    def calc_avg_bond_length(self, bond_lengths):
        sum_bond = 0
        counter = 0
        for i in range(0, len(bond_lengths)):
            for j in range(0, len(bond_lengths[i])):
                sum_bond += float(bond_lengths[i][j])
                counter += 1
        return sum_bond/counter

    def get_effective_cn(self, bond_lengths):
        ecn = []
        for i in bond_lengths:
            i.sort()
            lmin = i[0]
            num = 0.0
            den = 0.0
            for j in i:
                exp = (math.exp(1-(j/lmin)**6))
                num += j*exp
                den += exp
            lav = num/den
            weight = 0.0
            for j in i:
                weight += math.exp(1-(j/lav)**6)
            ecn.append(weight)
        return ecn

    def get_distortion_index(self, bond_lengths):
        di = []
        for i in bond_lengths:
            i.sort()
            lmin = i[0]
            num = 0.0
            den = 0.0
            for j in i:
                exp = (math.exp(1-(j/lmin)**6))
                num += j*exp
                den += exp
            lav = num/den
            DI = 0.0
            for j in i:
                DI += abs(j-lav)/lav
            di.append(DI/len(i))
        return di

    def get_path_length(self, structures, diffuser_index):
        # This should correctly handle any periodic boundary crossings along the diffusion pathway, but use caution.
        path_distance = 0
        for i in range(0, len(structures)-1):
            path_distance = path_distance + structures[i][diffuser_index].distance(structures[i+1][diffuser_index])
        return path_distance

    def compose_vasp_movie(self, diffusion_sites, empty_lattice):
        for i in diffusion_sites:
            site_from_dict = Site.from_dict(i)
            species = str(site_from_dict.specie)
            coords = site_from_dict.coords
            empty_lattice.append(species, coords, coords_are_cartesian=True)
        # writer.write_structure(empty_lattice, "POSCAR.path.vasp")
        return

    def render_path(self):
        directories = self.dirlist(os.getcwd())
        dirs_num = self.sort_dirs(directories)
        structure_array = self.vasp_neb_to_pymat(dirs_num)
        host_structure = self.get_empty_lattice(structure_array, self._working_ion)
        diffuser_index = self.get_diffuser_index(structure_array, self._working_ion)
        diffusion_sites = self.get_path_sites(structure_array, diffuser_index)
        self.compose_vasp_movie(diffusion_sites, host_structure)
        return

    def collect_inputs(self, root_path="./"):
        incar = Incar.from_file(root_path + "/INCAR")
        kpoints = Kpoints.from_file(root_path + "/KPOINTS")
        potcar = Potcar.from_file(root_path + "/POTCAR")
        potcar_array = {}
        for atom in potcar:
            potcar_array[atom.keywords["TITEL"]] = atom.keywords["VRHFIN"]
        input_dict = {"INCAR": incar.as_dict(), "KPOINTS": kpoints.as_dict(), "Pseudos": potcar_array}
        return input_dict

    def check_volume(self, structure_array):
        # Do we need to define a tolerance? Why not, this will probably save us time in the long run.
        tolerance = 0.001
        return (structure_array[0].volume - structure_array[-1].volume) < tolerance

    def populate_fields(self):
        directories = self.dirlist(os.getcwd())
        dirs_num = self.sort_dirs(directories)
        path_id = self._mp_id + "-NEB-" + self._pair_index

        structure_array = self.vasp_neb_to_pymat(dirs_num)
        structure_db = []
        for struct in structure_array:
            structure_db.append(struct.as_dict())
        poscar_array = self.get_all_poscars(dirs_num)
        host_structure = self.get_empty_lattice(structure_array, self._working_ion)

        all_cart_list = self.get_cart_images(dirs_num)
        atommap = self.get_atommap(dirs_num, poscar_array)

        cation_concentration = self.get_cation_concentration(structure_array, self._working_ion)

        diffuser_index = self.get_diffuser_index(structure_array, self._working_ion)
        diffusing_species = self.get_actual_diff_sp(poscar_array, all_cart_list, atommap, user_sp=self._working_ion,
                                                    user_pos=(diffuser_index+1))

        cation_diffusion_start, cation_diffusion_end = self.get_start_end_sites(structure_array, diffuser_index)
        diffusion_sites = self.get_path_sites(structure_array, diffuser_index)
        length = self.get_path_length(structure_array, diffuser_index)

        cn_path = self.get_cn(poscar_array, diffusing_species)
        voronoi_sites = self.get_voronoi(poscar_array, diffusing_species)
        neighbor_array = self.get_neighbors(poscar_array, atommap, voronoi_sites)
        bond_length_array = self.get_nn_bond_length(poscar_array, neighbor_array, diffusing_species)
        image_energies = self.get_energies_images(dirs_num)
        path_energies = self.get_scaled_energies(image_energies)
        activation_energy = self.max_energy(path_energies)
        avg_bond_length = self.calc_avg_bond_length(bond_length_array)
        ecn_path = self.get_effective_cn(bond_length_array)
        vol_relax_flag = str(self.check_volume(structure_array))
        input_set = self.collect_inputs(os.getcwd())

        distortion_index = self.get_distortion_index(bond_length_array)

        # Create dictionary for database fields
        database_fields = {'path_id': path_id, 'cation_concentration': cation_concentration, 'path_type': "NEB",
                           'host_structure': host_structure.as_dict(),
                           'cation_diffusion_start': cation_diffusion_start.as_dict(),
                           'cation_diffusion_end': cation_diffusion_end.as_dict(),
                           'path_sites': diffusion_sites, 'path_energy': path_energies.tolist(), 'path_length': length,
                           'coordination': cn_path, 'effective_coordination': ecn_path,
                           'average_bond_length': avg_bond_length, 'distortion_index': distortion_index,
                           'activation_energy': activation_energy, 'neb_images': structure_db,
                           'vol_relaxed': vol_relax_flag, 'input_set': input_set}
        return database_fields

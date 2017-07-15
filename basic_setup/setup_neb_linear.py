#!/usr/bin/python
try:
    import sys
    from pymatgen import Structure, Lattice
    from pymatgen.io.vasp           import Poscar
    from pymatgen.alchemy.transmuters import PoscarTransmuter
    from pymatgen.alchemy.transmuters import CifTransmuter
    from pymatgen.io.vasp.sets import MITRelaxSet 
except ImportError:
    print ("Unable to import pymatgen module.")
    sys.exit()

import os
import sys
import argparse
import functools
import numpy as np



def makefolders(root_dir, subfolders):
    concat_path = functools.partial(os.path.join, root_dir)
    map(os.makedirs, map(concat_path, subfolders))

def interpolate(s1, s2, nimages):
    '''
    doesn't relax either end
    '''
    fcoords1 = np.array(s1.frac_coords)
    fcoords2 = np.array(s2.frac_coords)
    lvect = s2.lattice.matrix - s1.lattice.matrix
    fcvect = fcoords2 - fcoords1 - np.round(fcoords2 - fcoords1)
    structures = []
    for x in xrange(nimages):
        l = Lattice(s1.lattice.matrix + x * lvect / (nimages - 1))
        fc = fcoords1 + x * fcvect / (nimages - 1)
        s = Structure(l, s1.species, fc)
        structures.append(s)
    return structures

def prepareFiles (sp, images):


    default_texas = { 
                     'IOPT'  : '1', 'IBRION' : '3', 
                     'POTIM' : '0', 'MAXMOVE' : '0.2' , 
                     'ILBFGSMEM' : '20', 'LGLOBAL' : '.TRUE.', 
                     'LAUTOSCALE' : '.TRUE.', 'INVCURV' : '0.01', 
                     'LLINEOPT' : '.FALSE.', 'FDSTEP' : '5E-3',
                     'LCLIMB' : '.FALSE.', 'SPRING' : '-5.0', 
                     'NSIM' : '4'  , 'NPAR' : '4', 
                     'LPLANE' : '.TRUE.' , 'LWAVE' : '.FALSE.', 
                     'LCHARG' : '.FALSE.', 'ICHARG' : '2', 
                     'LREAL' : '.FALSE.', 'ENCUT' : '500',
                     'ISIF'  : 2, 'EDIFF': 1e-5, 'LDAU': '.FALSE.', 
                     'NSW' : 200, 'IMAGES' : images, 'EDIFFG' : -0.05
                     } 
   
    vasp_set = MITRelaxSet(sp, user_incar_settings=default_texas)
    vasp_set.write_input("./neb/")  

    

    
parser = argparse.ArgumentParser()

parser.add_argument(
        "images",
        help="No. images")


args = parser.parse_args()

images = int(args.images)

'''
intial and last images
'''
sp = Structure.from_file("./endpt_relax/initial/CONTCAR")
ep = Structure.from_file("./endpt_relax/final/CONTCAR")

os.makedirs("./neb/")

prepareFiles(sp, images)


for i, s in enumerate(interpolate(sp, ep, images + 2)):
    os.makedirs("./neb/" + os.path.join(str(i).zfill(2)))
    s.to(fmt="poscar", filename=("./neb/" + str(i).zfill(2) + "/POSCAR"))
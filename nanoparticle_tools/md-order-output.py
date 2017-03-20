import sys
import csv
import math
import scipy.special
import scipy.misc
import numpy as np

# Need me some subroutines! Mmmm, subroutines.  atomData[4][i][1] should be the polar angle, theta, 
# and atomData[4][i][2] should be the azimuthal angle, phi.  SciPy reverses these.
def computeQLM(m,l,atomData):
  sum = 0
  #print "m is " + str(m)
  #print "l is " + str(l)
  for i in range (0, len(atomData[4])):
    #print "Polar angle for neighbor " + str(i) + " is " + str(atomData[4][i][1])
    #print "Azimuthal angle for neighbor " + str(i) + " is " + str(atomData[4][i][2])
    sum = sum + scipy.special.sph_harm(m,l,atomData[4][i][2],atomData[4][i][1])
  if len(atomData[4]) > 0:
    sum = sum/(len(atomData[4]))
  #print sum
  return sum
  

# Start by opening a .csv with a given separator - for LAMMPS, this is a space.

# Process command-line arguments: 
#	1.  Filename (XYZ file with first two lines removed)
#	2.  Bond length
#	3.  Element in question (case-sensitive, probably)
#	4.  The q order to compute
#	5.  Output file
xyz_file = sys.argv[1]
cutoff = float(sys.argv[2])
element = sys.argv[3]
qOrder = int(sys.argv[4])
out_file = open(sys.argv[5], 'w')

# Open the xyz file, create a list of 4-element arrays structured in the following way:
#	atomList = [atom1, atom2, ..., atomN]
#	atom1 = ["Ge", x-coord, y-coord, z-coord]
#	atom2 = ["Ge", x-coord, y-coord, z-coord]
# And so forth.
separator = ' '
coordReader = csv.reader(open(xyz_file, 'rb'), delimiter=separator)
coordReader.next()
coordReader.next()
atomList = []

for row in coordReader:
  atomName = row[0]
  xcoord = row[1]
  ycoord = row[2]
  zcoord = row[3]
  atomData = []
  if atomName == element:
    atomData.append(atomName)
    atomData.append(float(xcoord))
    atomData.append(float(ycoord))
    atomData.append(float(zcoord))
    atomList.append(atomData)
  
# Build a nearest neighbor list.  These are stored by appending r, theta, phi to the atomData for that atom:
#	atomList[N] = ["Ge", xcoord, ycoord, zcoord, [[r1, theta1, phi1, j1], [r2, theta2, phi2, j2], ...,[rN, thetaN, phiN, jN]]
# And so forth.  The number of nearest neighbors for atom i is given by len(atomList[i][4]).
for i in range (0, len(atomList)):
  neighbors = []
  for j in range (0, len(atomList)):
    if atomList[i] != atomList[j]:
      dx = atomList[i][1] - atomList[j][1]
      dy = atomList[i][2] - atomList[j][2]
      dz = atomList[i][3] - atomList[j][3]
      rsq = dx*dx + dy*dy + dz*dz
      r = math.sqrt(rsq)
      bondData = []
      if r <= cutoff:
	theta = math.acos(dz/r)
	phi = math.acos(dx/r*math.sin(theta))
	bondData.append(r)
	bondData.append(theta)
	bondData.append(phi)
	bondData.append(j)
        neighbors.append(bondData)
  atomList[i].append(neighbors)
  #print "Neighbor vector for atom i is " + str(atomList[i][4])
 
# Compute the spherical harmonic vectors for atoms and their nearest neighbors.
# The efficiency of this program could be dramatically improved by combining these loops.
# For clarity, I'm leaving them separate for the moment.
for i in range (0, len(atomList)):
  qlVector = []
  for j in range (0, 2*qOrder):
    qlmSum = computeQLM(-qOrder+j, qOrder, atomList[i])
    qlVector.append(qlmSum)
  atomList[i].append(qlVector)
  
# At this point in the code, every element of atomList has a qlVector added to it as element [5].
# The accessors are about to get crazy when we dot these vectors together.
# Let's compose q for each atom.  The index of neighbor j for atom i is atomList[i][4][j][3].
# Therefore the QL vector for the neighbor is given by atomList[atomList[i][4][j][3]][5].
histoListo = []
for i in range (0, len(atomList)):
  sum = 0.00 + 0.00j
  for j in range (0, len(atomList[i][4])):
    sum = sum + np.vdot(atomList[atomList[i][4][j][3]][5], atomList[i][5])/(np.linalg.norm(atomList[i][5])*np.linalg.norm(np.ma.conjugate(atomList[atomList[i][4][j][3]][5])))
  if len(atomList[i][4]) > 0:
    sum = sum/len(atomList[i][4])
  histoListo.append(sum)
  
# Generate output files
for i in range(0, len(atomList)):
    out_file.write(str(atomList[i][0]) + "\t" + str(len(atomList[i][4])) + "\t" + str(histoListo[i].real) + "\n")

import sys
import csv
import math
import os
import numpy as np

def getNumbers(trajectory_file):
	atomsPerFrame = 0
	linesInFile = sum(1 for line in open(trajectory_file))
	with open(trajectory_file) as file:
		atomsPerFrame = int(file.readline())
	file.close()
	return atomsPerFrame, (linesInFile/(atomsPerFrame + 2))

# This routine is for hexylamine passivated CdSe and is designed to work on one frame of an XYZ trajectory.
def countMolecules(trajectory_frame, totalAtoms, atomsPerSolvent, atomsPerLigand):
	numCd = 0
	numSe = 0
	numLigand = 0
	for i in range(0, len(trajectory_frame)):
		numCd = numCd + trajectory_frame[i].count('Cd')
		numSe = numSe + trajectory_frame[i].count('Se')
		numLigand = numLigand + trajectory_frame[i].count('HN')
	numLigand = numLigand/2
	numSolvent = ((totalAtoms) - (numCd + numSe + numLigand*atomsPerLigand))/atomsPerSolvent
	return numCd, numSe, numLigand, numSolvent

def py_angle(vector1, vector2):
	# We use the tan function in this case because cos yields inaccurate results for small angles.
	cosang = np.dot(vector1, vector2)
	sinang = np.linalg.norm(np.cross(vector1, vector2))
	return np.arctan2(sinang, cosang)

def computeQab(eiAarray, eiBarray):
	deltaAB = 0
	if eiAarray == eiBarray:
		deltaAB = 1
	sum = 0
	for i in range(0, len(eiAarray)):
		sum = sum + (3*eiAarray[i]*eiBarray[i] - deltaAB)
	avg = sum/(2.0*len(eiAarray))
	return avg

def computeEI(split_array):
	ei = []
	eix = float(split_array[0][1]) - float(split_array[len(split_array)-1][1])
	eiy = float(split_array[0][2]) - float(split_array[len(split_array)-1][2])
	eiz = float(split_array[0][3]) - float(split_array[len(split_array)-1][3])
	ei.append(eix)
	ei.append(eiy)
	ei.append(eiz)
	ei_magnitude = math.sqrt(ei[0]*ei[0] + ei[1]*ei[1] + ei[2]*ei[2])
	ei_unit = []
	ei_unit.append(float(ei[0]/ei_magnitude))
	ei_unit.append(float(ei[1]/ei_magnitude))
	ei_unit.append(float(ei[2]/ei_magnitude))
	return ei_unit

def makevectorarray(eiList):
	xArray = []
	yArray = []
	zArray = []
	for i in range(0, len(eiList)):
		xArray.append(eiList[i][0])
		yArray.append(eiList[i][1])
		zArray.append(eiList[i][2])
	return xArray, yArray, zArray

def buildQtensor(eiArray):
	eixArray, eiyArray, eizArray = makevectorarray(eiArray)
	Qxx = computeQab(eixArray, eixArray)
	Qxy = computeQab(eixArray, eiyArray)
	Qxz = computeQab(eixArray, eizArray)
	Qyx = computeQab(eiyArray, eixArray)
	Qyy = computeQab(eiyArray, eiyArray)
	Qyz = computeQab(eiyArray, eizArray)
	Qzx = computeQab(eizArray, eixArray)
	Qzy = computeQab(eizArray, eiyArray)
	Qzz = computeQab(eizArray, eizArray)
	Qtensor = np.array([[Qxx, Qxy, Qxz], [Qyx, Qyy, Qyz,], [Qzx, Qzy, Qzz]])
	return Qtensor

# Getting the file into an array of frames
numAtoms, numFrames = getNumbers(sys.argv[1])
separator = ' '
coordReader = csv.reader(open(sys.argv[1], 'rb'), delimiter=separator)
frameArray = []
i = 0
j = 0
while i < numFrames:
	frameData = []
	coordReader.next()
	coordReader.next()
	while j < (numAtoms):
		rowHolder = coordReader.next()
		atomName = rowHolder[0]
		xcoord = rowHolder[1]
		ycoord = rowHolder[2]
		zcoord = rowHolder[3]
		atomData = []
		atomData.append(atomName)
		atomData.append('%f' % float(xcoord))
		atomData.append('%f' % float(ycoord))
		atomData.append('%f' % float(zcoord))
		frameData.append(atomData)
		j = j + 1
	frameArray.append(frameData)
	j = 0
	i = i + 1

# We still need some numbers to work with
atomsPerSolvent = 6
atomsPerLigand = 9
numCd, numSe, numLigand, numSolvent = countMolecules(frameArray[0], numAtoms, atomsPerSolvent, atomsPerLigand)

# Now, for each frame, we want to compute the director axis.
dotArray = []
for i in range(0, len(frameArray)):
	cdse_array = []
	ligand_array = []
	solvent_array = [] 
	for j in range(0, len(frameArray[i])):
		print "Analyzing frame " + str(i)
		if j < (numCd + numSe):
			cdse_array.append(frameArray[i][j])
		elif j < (numCd + numSe + numLigand*atomsPerLigand):
			ligand_array.append(frameArray[i][j])
		else:
			solvent_array.append(frameArray[i][j])
	ligand_molecule_array = np.split(np.array(ligand_array), numLigand)
	solvent_molecule_array = np.split(np.array(solvent_array), numSolvent)
	ligand_ei_array = []
	solvent_ei_array = []
	for k in range(0, len(ligand_molecule_array)):
		ei = computeEI(ligand_molecule_array[k])
		ligand_ei_array.append(ei)
	for k in range(0, len(solvent_molecule_array)):
		ei = computeEI(solvent_molecule_array[k])
		solvent_ei_array.append(ei)
	ligandQtensor = buildQtensor(ligand_ei_array)
	solventQtensor = buildQtensor(solvent_ei_array)
	ligandQeigval, ligandQeigvec = np.linalg.eig(ligandQtensor)
	solventQeigval, solventQeigvec = np.linalg.eig(solventQtensor)
	ligandDirectorAxis = ligandQeigvec[:,ligandQeigval.argmax(axis=0)]
	solventDirectorAxis = solventQeigvec[:,solventQeigval.argmax(axis=0)]
	directorAngle = 57.2957795*py_angle(ligandDirectorAxis, solventDirectorAxis)
	if directorAngle > 90.0:
		ligandDirectorAxis = -1.0*ligandDirectorAxis
		directorAngle = 57.2957795*py_angle(ligandDirectorAxis, solventDirectorAxis)
	directorDot = np.dot(ligandDirectorAxis, solventDirectorAxis)
	xy_pair = [i, directorDot]
	dotArray.append(xy_pair)

outfile = open(sys.argv[2], 'w')
for i in range(0, len(dotArray)):
	outfile.write(str(dotArray[i][0]) + "\t" + str(dotArray[i][1]) + "\n")


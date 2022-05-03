from .molecule import Molecule
import sys
write = sys.stdout.write

def read_pdb(filename):

	atomsym = []
	xyz = []
	pdbf = open(filename,"r")
	for line in pdbf.readlines():
		columns = line.split()
		ncols = len(columns)
		if ncols >= 11 and columns[0] == 'ATOM':
			atomsym.append(columns[2][0:1])
			xyz.append(float(columns[-7]))
			xyz.append(float(columns[-6]))
			xyz.append(float(columns[-5]))

	molecule = Molecule()
	molecule.setup(atomsym, xyz)

	return molecule


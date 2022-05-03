from .molecule import Molecule
from .util import is_integer, atmsym
import sys
write = sys.stdout.write

def read_xyz(filename):

	atomsym = []
	xyz = []
	xyzf = open(filename,"r")
	for line in xyzf.readlines():
		columns = line.split()
		ncols = len(columns)
		if ncols == 4:
			atomsym.append(columns[0])
			xyz.append(float(columns[1]))
			xyz.append(float(columns[2]))
			xyz.append(float(columns[3]))

	molecule = Molecule()
	molecule.setup(atomsym, xyz)

	return molecule

def read_one_xyz_out_of_many(filename, imolecule):

	atomsym = []
	xyz = []
	xyzf = open(filename, "r")
	#figure out number of atoms and gap lines
	natoms = 0
	ngap   = 0
	line1  = 0
	icount = 0
	for line in xyzf.readlines():
		columns = line.split()
		ncols = len(columns)
		if ncols == 1: 
			if natoms == 0: 
				natoms = int(columns[0])
				line1  = icount
			elif is_integer(columns[0]):
				ngap   = icount-line1-natoms-2
				break
		icount += 1
	#print("natoms: ", natoms, "ngap:", ngap)
	xyzf.close()
	xyzf = open(filename, "r")
	i = -1
	for line in xyzf.readlines():
		columns = line.split()
		ncols = len(columns)
		if ncols == 1 and i == -1: i = 1
#		if ncols == 1 and natoms == 0: 
#			natoms = int(columns[0])
#			#write("read in coordinates of %3d atoms\n" % natoms)
#			i = 1
#		elif i >= (imolecule-1)*(natoms+2) + 2 and i <= imolecule*(natoms+2) and ncols == 4:
		elif i >= (imolecule-1)*(natoms+2+ngap) + 2 and i <= imolecule*(natoms+2+ngap) and ncols == 4:
			if columns[0] == 'OH2':
				atomsym.append('O')
			elif columns[0] == 'H1' or columns[0] == 'H2':
				atomsym.append('H')
			elif is_integer(columns[0]): 
				atomsym.append(atmsym(int(columns[0])))
			else:
				atomsym.append(columns[0])
			xyz.append(float(columns[1]))
			xyz.append(float(columns[2]))
			xyz.append(float(columns[3]))

		i += 1

	molecule = Molecule()
	molecule.setup(atomsym, xyz)

	return molecule

def read_all_geometries(xyzfile, natoms):

    xyzf = open(xyzfile, "r")
    nlines = len(xyzf.readlines())
    nframes = int(nlines/(natoms+2))
    write("Reading %d frames from the IRC file\n" %(nframes))
    frames = []
    for iframe in range(0, nframes):
            frame = read_one_xyz_out_of_many(xyzfile, iframe+1)
            frames.append(frame)
    xyzf.close()
    return frames

def compare_molecules(mol1, mol2):
        tol = 0.0001
        is_same = True
        natoms1 = mol1.natoms 
        natoms2 = mol2.natoms 
        print("natoms:", natoms1, natoms2)
        if natoms1 != natoms2: is_same = False
        for k in range(0, natoms1):
                for m in range(0, 3): 
                        if abs(mol1.coords[3*k+m] - mol2.coords[3*k+m]) > tol:
                                #print(k, m, mol1.coords[3*k+m], mol2.coords[3*k+m])
                                is_same = False
        return is_same


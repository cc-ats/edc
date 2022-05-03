import sys 
import numpy as np
write = sys.stdout.write
from .cubic_grid import *
from .util import *

def write_cube(filename, molecule, grid, grid_data):

	write("writing to %s\n" % filename)
	cubef = open(filename, "w")
	cubef.write("filename: %s\n\n" % (filename))
	cubef.write("%5d %12.6f %12.6f %12.6f\n" % (molecule.natoms, grid.x[0], grid.y[0], grid.z[0]))
	cubef.write("%5d %12.6f %12.6f %12.6f\n" % (grid.nx, grid.xinc, 0.0, 0.0))
	cubef.write("%5d %12.6f %12.6f %12.6f\n" % (grid.ny, 0.0, grid.yinc, 0.0))
	cubef.write("%5d %12.6f %12.6f %12.6f\n" % (grid.nz, 0.0, 0.0, grid.zinc))
	for k in range(0, molecule.natoms):
		cubef.write("%5d %12.6f" % (molecule.atmnum[k], molecule.atmnum[k]*1.0))
		for m in range(0, 3): 
			cubef.write(" %12.6f" % (molecule.coords[3*k+m]))
		cubef.write("\n")
	count = 0 
	for ix in range(0, grid.nx):
		for iy in range(0, grid.ny):
			for iz in range(0, grid.nz):
				cubef.write(" %16.9E" % grid_data[count])
				count += 1
				if int((iz+1)/6)*6 == iz+1: cubef.write("\n")
			if int(grid.nz/6)*6 != grid.nz: cubef.write("\n");
	cubef.close()

def write_cube_squared(filename, molecule, grid, grid_data):

	write("writing to %s\n" % filename)
	cubef = open(filename, "w")
	cubef.write("filename: %s\n\n" % (filename))
	cubef.write("%5d %12.6f %12.6f %12.6f\n" % (molecule.natoms, grid.x[0], grid.y[0], grid.z[0]))
	cubef.write("%5d %12.6f %12.6f %12.6f\n" % (grid.nx, grid.xinc, 0.0, 0.0))
	cubef.write("%5d %12.6f %12.6f %12.6f\n" % (grid.ny, 0.0, grid.yinc, 0.0))
	cubef.write("%5d %12.6f %12.6f %12.6f\n" % (grid.nz, 0.0, 0.0, grid.zinc))
	for k in range(0, molecule.natoms):
		cubef.write("%5d %12.6f" % (molecule.atmnum[k], molecule.atmnum[k]*1.0))
		for m in range(0, 3): 
			cubef.write(" %12.6f" % (molecule.coords[3*k+m]))
		cubef.write("\n")
	count = 0 
	for ix in range(0, grid.nx):
		for iy in range(0, grid.ny):
			for iz in range(0, grid.nz):
				cubef.write(" %16.9E" % (grid_data[count]*grid_data[count]))
				count += 1
				if int((iz+1)/6)*6 == iz+1: cubef.write("\n")
			if int(grid.nz/6)*6 != grid.nz: cubef.write("\n");
	cubef.close()

def read_cube_file(filename):
	cubef = open(filename, "r")
	lines = cubef.readlines()
	natoms = int(lines[2].split()[0])
	print("natoms: ", natoms)

	#grid information
	ni = []
	cmin = []
	cinc = []
	ni.append(int(lines[3].split()[0]))
	ni.append(int(lines[4].split()[0]))
	ni.append(int(lines[5].split()[0]))
	cmin.append(float(lines[2].split()[1]))
	cmin.append(float(lines[2].split()[2]))
	cmin.append(float(lines[2].split()[3]))
	cinc.append(float(lines[3].split()[1]))
	cinc.append(float(lines[4].split()[2]))
	cinc.append(float(lines[5].split()[3]))
	grid = CubicGrid()
	grid.setup2(ni, cmin, cinc)

	#molecule information
	atmnum = []
	coord = []
	for k in range(0, natoms):
		cols = lines[6+k].split()
		atmnum.append(int(cols[0]))
		coord.append(float(cols[2])*BOHR)
		coord.append(float(cols[3])*BOHR)
		coord.append(float(cols[4])*BOHR)
	mol = Molecule()
	mol.setup2(atmnum, coord)
	mol.mprint("This is the molecule")

	#actual data
	grid_data = []
	for k in range(natoms+6, len(lines)):
		cols = lines[k].split()
		ncols = len(cols)
		for m in range(0, ncols):
			grid_data.append(float(cols[m]))

	return mol, grid, grid_data

def scale_cube_data(filename, filename2, scale):
	#read information
	mol, grid, grid_data = read_cube_file(filename)

	#scale grid_data
	grid_data = np.multiply(scale,  grid_data)

	#write to cube file
	write_cube(filename2, mol, grid, grid_data)
	

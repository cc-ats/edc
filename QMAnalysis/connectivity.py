import numpy as np
from math import *
from .util import *
from .molecule import *
import sys
write = sys.stdout.write

def print_distance_and_connectivity(distances, connectivity, molecule):
	write("Covalent bonds in the molecule:\n")
	natoms = molecule.natoms
	for k in range(0, natoms):
		for l in range(k+1, natoms):
			if connectivity[k+l*natoms] == 1:
				write("   %2s%-3d - %2s%-3d: %9.5f Ang\n" % 
				      (molecule.atmsym[k], k+1, molecule.atmsym[l], l+1, distances[k+l*natoms]))
	write("\n")

def print_distances_2_img(xyzfile, distances, connectivity, molecule):

	spt = open("connect.spt", "w")
	spt.write("background = white\n")
	spt.write("load "+xyzfile+"\n")
	spt.write("label %a\n")
	spt.write("defaultDistanceLabel = \"%6.4VALUE\"\n")
	spt.write("set measurements angstroms\n")
	spt.write("set measurementnumbers on\n")
	natoms = molecule.natoms
	for k in range(0, natoms):
		for l in range(k+1, natoms):
			if connectivity[k+l*natoms] == 1 and molecule.atmsym[k] != 'H' and molecule.atmsym[l] != 'H':
				spt.write("measure (atomno=%3d) (atomno=%3d)\n" % (k+1, l+1))
	spt.write("write image "+xyzfile[:-3]+"png\n")
	spt.close()
	os.system("jmol -n connect.spt")

def get_distance_and_connectivity(distances, connectivity, molecule):

	natoms = molecule.natoms
	#write("natoms=%d\n" % (natoms))

	connect = []
	for k in range(0, natoms*natoms):
		connect.append(0)

	# nearest neighbors
	for i in range(0, natoms):
		for j in range(i+1, natoms):
			#compute distance
			r2 = 0
			for k in range(0, 3):
				diff = molecule.coords[3*i+k] - molecule.coords[3*j+k]
				r2 += diff*diff
			dist = sqrt(r2)*BOHR
			distances[i+j*natoms] = distances[j+i*natoms] = dist

			#connectivitity
			#this is a little bit arbitrary, need to fix later
			is_bond = 0
			if molecule.atmsym[i] != 'H' and molecule.atmsym[j] != 'H':
				is_bond = dist < 1.9
			elif molecule.atmsym[i] != 'H' or molecule.atmsym[j] != 'H':
				if molecule.atmsym[i] == 'O' or molecule.atmsym[j] == 'O':  #O-H bond
					is_bond = dist < 1.2
				else:
					is_bond = dist < 1.6
			else:
				is_bond = dist < 1.2
			if is_bond: 
				connect[i+j*natoms] = connect[j+i*natoms] = 1

        # second-nearest neighbors
	for i in range(0, natoms):
		for j in range(0, natoms):
			if connect[i+j*natoms] == 1:
				for k in range(0, natoms):
					if k != i and k != j and connect[k+j*natoms] == 1 and connect[k+i*natoms] == 0:
						connect[k+i*natoms] = connect[i+k*natoms] = 2

	# third-nearest neighbors
	for i in range(0, natoms):
		for j in range(0, natoms):
			if connect[i+j*natoms] == 2:
				for k in range(0, natoms):
					if k != i and k != j and connect[k+j*natoms] == 1 and connect[k+i*natoms] == 0:
						connect[k+i*natoms] = connect[i+k*natoms] = 3

	connectivity[:] = connect

def fragment_analysis(frag_index, natoms, connectivity):
	ifrag = -1
	for m in range(0, natoms): frag_index.append(-1)
	for m in range(0, natoms): 
		iatom = -1
		for n in range(0, natoms):
			if frag_index[n] == -1:	
				iatom = n	
				ifrag += 1
				break
		#print("iatom: ", iatom, "ifrag=", ifrag)
		if iatom == -1: break
		frag_index[iatom] = ifrag
		for n in range(0, natoms):
			if frag_index[n] == -1:
				for k in range(0, natoms):
					if frag_index[k] == ifrag and connectivity[k+n*natoms] == 1:
						frag_index[n] = ifrag
	return ifrag+1



# -*- coding: utf-8 -*-
#
# Copyright (c), the QMAnalysis development team
#
#
#

""" Parser for Q-Chem output file """

import numpy, sys, os
from .connectivity import *
from .util import *
from .math_util import *
write = sys.stdout.write

class Molecule():
	""" Q-Chem output file. """

	def __init__(self):
		#super(QChem, self).__init__(logname="QChem")
		#doing nothing
		if 1 == 0: write("initialization")

	def setup(self, atmsym, coords):
		"""setup molecule with atomic symbols"""
		self.atmsym = atmsym
		self.coords = coords
		self.natoms = len(atmsym)
		self.atmnum = []
		for k in range(0, self.natoms):
			self.coords[3*k+0] *= 1/BOHR
			self.coords[3*k+1] *= 1/BOHR
			self.coords[3*k+2] *= 1/BOHR
			self.atmnum.append(atmnum(atmsym[k]))
		self.charge = 0
		self.multiplicity = 1

	def setup2(self, atmnum, coords):
		"setup molecule with atomic numbers"
		self.atmnum = atmnum
		self.coords = coords
		self.natoms = len(atmnum)
		self.atmsym = []
		for k in range(0, self.natoms):
			self.coords[3*k+0] *= 1/BOHR
			self.coords[3*k+1] *= 1/BOHR
			self.coords[3*k+2] *= 1/BOHR
			self.atmsym.append(atmsym(atmnum[k]))
		self.charge = 0
		self.multiplicity = 1

	def setup_submolecule(self, molecule2, atom_list):
		"""setup molecule using a subset of atoms from another molecule"""
		self.atmsym = []
		self.atmnum = []
		self.coords = []
		for m in range(0, len(atom_list)):
			k = atom_list[m] - 1
			self.atmsym.append(molecule2.atmsym[k])
			self.atmnum.append(molecule2.atmnum[k])
			self.coords.append(molecule2.coords[3*k+0])
			self.coords.append(molecule2.coords[3*k+1])
			self.coords.append(molecule2.coords[3*k+2])
		self.natoms = len(self.atmsym)
		self.charge = 0
		self.multiplicity = 1

	def setup_z(self, z_elements, z_values):
		self.natoms = len(z_elements)
		self.z_elements = z_elements
		self.z_values   = z_values
		self.charge = 0
		self.multiplicity = 1

	def set_charge_multiplicity(self, charge, multiplicity):
		self.charge = charge
		self.multiplicity = multiplicity

	def mprint(self, title):
		write("  ===== %s (in Angstrom) =====\n" % (title))
		for k in range(0, self.natoms):
			write ("%3s " % (self.atmsym[k]))
			for m in range(0, 3):
				write ("%12.7f " % (self.coords[3*k+m]*BOHR))
			write ("\n")
		write ("\n")
	
	def mprint_with_energy(self, xyzfile, energy, append='none'):
		if append == 'none':
			xyzf = open(xyzfile, "w")
		else:
			xyzf = open(xyzfile, "a")
		xyzf.write("%d \n" % self.natoms)
		xyzf.write("%15.10f \n" % energy)
		for k in range(0, self.natoms):
			xyzf.write ("%3s " % (self.atmsym[k]))
			for m in range(0, 3):
				xyzf.write ("%12.7f " % (self.coords[3*k+m]*BOHR))
			xyzf.write("\n")

	def molecular_fragment_analysis(self):
		#get connectivity
		self.distances = np.zeros(self.natoms*self.natoms)
		self.connectivity = []
		for k in range(0, self.natoms*self.natoms): self.connectivity.append(0)
		get_distance_and_connectivity(self.distances, self.connectivity, self)
		#print("connectivity=", connectivity)

		#fragment analysis
		self.frag_index = []
		self.nfrag = fragment_analysis(self.frag_index, self.natoms, self.connectivity)
		#write("nfrag= %d\n" % self.nfrag)
		#for k in range(0, self.natoms):
		#	write("frag_index= %d\n" % self.frag_index[k])

	def find_n_closest_fragments(self, nmax):
		#relative to first fragment
		nfrag = self.nfrag
		fragment_distances = np.zeros(nfrag)
		for i in range(0, nfrag): fragment_distances[i] = 1000.0
		for i in range(0, self.natoms):
			if self.frag_index[i] == 0:
				for j in range(0, self.natoms):
					jfrag = self.frag_index[j] 
					if jfrag != 0:
						if self.distances[i+j*self.natoms] < fragment_distances[jfrag]:
							fragment_distances[jfrag] = self.distances[i+j*self.natoms]
		matrix_print_1d(fragment_distances, 1, nfrag, 6, "fragment distances") 
		indices = np.argsort(fragment_distances)
		print("indices=", indices)
		atomlist = []
		for i in range(0, self.natoms):
			ifrag = self.frag_index[i]
			if (ifrag == 0):
				atomlist.append(i+1)
			else:
				for m in range(0, nmax):
					if indices[m] == ifrag:
						atomlist.append(i+1)
		print("atomlist: ", atomlist)
		return atomlist

	def get_connectivity(self):
		#get connectivity
		self.distances = np.zeros(self.natoms*self.natoms)
		self.connectivity = []
		for k in range(0, self.natoms*self.natoms): self.connectivity.append(0)
		get_distance_and_connectivity(self.distances, self.connectivity, self)
		#print("connectivity=", connectivity)

	def stretch_one_bond(self, atom1, atom2, dist):
		self.get_connectivity()
		natoms = self.natoms
		#matrix_int_print_2d(self.connectivity, natoms, natoms, "connectivity")
		connectivity2 = np.zeros(natoms*natoms)
		for k in range(0, natoms*natoms):
			connectivity2[k] = self.connectivity[k]
		connectivity2[atom1-1+(atom2-1)*natoms] = 0
		connectivity2[atom2-1+(atom1-1)*natoms] = 0
		print("atom1: ", atom1, "atom2:", atom2)
		#matrix_int_print_2d(connectivity2, natoms, natoms, "connectivity 2")
		frag_index = []
		nfrag = fragment_analysis(frag_index, natoms, connectivity2)
		print("nfrag=", nfrag)
		for k in range(0, natoms):
			print("k", k+1, "frag_index:", frag_index[k])
		index2 = frag_index[atom2-1]
		
		dist12 = self.distances[atom1-1+(atom2-1)*natoms]
		scale = dist/dist12 - 1.0
		print("dist12=", dist12, "scale=", scale)

		v12 = np.zeros(3)
		for k in range(0, 3):
			v12[k] = (self.coords[3*(atom2-1)+k] - self.coords[3*(atom1-1)+k])*scale

		coords = np.zeros(3*natoms)
		for i in range(0, natoms):
			for k in range(0, 3):
				coords[3*i+k] = self.coords[3*i+k] * BOHR
				if frag_index[i] == index2:
					coords[3*i+k] += v12[k] *BOHR

		molecule2 = Molecule()
		molecule2.setup(self.atmsym, coords)
		return molecule2



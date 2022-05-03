# -*- coding: utf-8 -*-
#
# Copyright (c), the QMAnalysis development team
#
#
#

import sys, os
import numpy as np
from .connectivity import *
from .molecule import *
from .util import *
write = sys.stdout.write

class ExternalCharges():

	def __init__(self):
		self.nchgs = 0
		if 1 == 0: write("initialization")

	def setup(self, coords, charges):
		self.nchgs   = len(charges)
		self.coords  = coords
		self.charges = charges
		self.cutoff_charges = np.zeros(self.nchgs)
		for k in range(0, self.nchgs):
			self.coords[3*k+0] *= 1/BOHR
			self.coords[3*k+1] *= 1/BOHR
			self.coords[3*k+2] *= 1/BOHR
			self.cutoff_charges[k] = self.charges[k]

	def reset_charges(self, new_charges):
		for k in range(0, self.nchgs):
			self.charges[k]        = new_charges[k]
			self.cutoff_charges[k] = new_charges[k]

	def copy(self, extCharges):
		self.nchgs = extCharges.nchgs
		self.coords  = extCharges.coords
		self.charges = extCharges.charges
		self.cutoff_charges = extCharges.cutoff_charges

	#it only works for water now
	def get_distances(self, molecule):
		self.distances = []
		for k in range(0, self.nchgs):
			if abs(self.charges[k] - 0.417) <0.001: continue    #water hydrogen
			nblock = 1
			if abs(self.charges[k] + 0.834) <0.001: nblock = 3  #water oxygen
			dmin = 10000.0
			for l in range(0, nblock):
				m = k + l
				for i in range(0, molecule.natoms):
					d2 = 0.0
					for j in range(0, 3):
						dj = self.coords[3*m+j] - molecule.coords[3*i+j] 
						d2 += dj * dj
					d = sqrt(d2)*BOHR
					if d < dmin: dmin = d
			for l in range(0, nblock):
				self.distances.append(dmin)


	#change charges with a cutoff criterion
	def set_cutoff_charges(self, cutoff):
		count = 0
		for k in range(0, self.nchgs):
			if self.distances[k] < cutoff:
				self.cutoff_charges[k] = self.charges[k]
				count += 1
			else: 
				self.cutoff_charges[k] = 0.0
		print(count, "charges within the cutoff distance of", cutoff, "Angstom")

	#atom site potential for charges beyonds the cutoff distance
	def compute_atom_site_potential(self, cutoff, molecule):
		natoms    = molecule.natoms
		potential = np.zeros(4*natoms)
		if not hasattr(self, 'distances'):
			self.get_distances(molecule)
		for k in range(0, self.nchgs):
			if self.distances[k] > cutoff:
				#print("k:", k, "charges=", self.charges[k])
				for i in range(0, molecule.natoms):
					dx = self.coords[3*k+0] - molecule.coords[3*i+0]
					dy = self.coords[3*k+1] - molecule.coords[3*i+1]
					dz = self.coords[3*k+2] - molecule.coords[3*i+2]
					d2 = dx*dx + dy*dy + dz*dz
					d  = sqrt(d2)
					potential[4*i+0] += self.charges[k] / d
					potential[4*i+1] += self.charges[k] * dx / (d*d2)
					potential[4*i+2] += self.charges[k] * dy / (d*d2)
					potential[4*i+3] += self.charges[k] * dz / (d*d2)
		return potential

				
	

# -*- coding: utf-8 -*-
#
# Copyright (c), the QMAnalysis development team
#
#
#


""" Parser for Q-Chem output file """

import sys, os
import numpy as np
from .util import *
from .molecule import Molecule
from .basis_set import *
from .tddft_functions import *

class QChemFChk():
	""" Q-Chem output file. """

	def __init__(self, fchkfile):
		#if 1 == 0: print "initialization of QChemFChk"
		self.fchkfile = fchkfile
		self.setup()

	def setup(self):
		self.natoms =int(get_first_value(self.fchkfile, "Number of atoms", -1))
		self.nalpha =int(get_first_value(self.fchkfile, "Number of alpha electrons", -1))
		self.nbeta  =int(get_first_value(self.fchkfile, "Number of beta electrons", -1))
		self.nbas   =int(get_first_value(self.fchkfile, "Number of basis functions", -1))
		self.norb   =int(int(get_first_value(self.fchkfile, "Alpha MO coefficients", -1))/self.nbas)
		self.noa    =self.nalpha
		self.nob    =self.nbeta
		self.nva    =self.norb-self.noa
		self.nvb    =self.norb-self.nob
		self.get_geometry()
		self.get_basis()

	def fchk_get_last_integer_values(self, phrase):
		values = []
		info = last_occurrence(self.fchkfile, phrase, 1)
		line1, nvalues = info[0], info[1]
		nlines = int(nvalues/6)
		if nlines * 6 < nvalues : nlines += 1
		line2 = line1 + nlines
		fchk = open(self.fchkfile, "r")
		for line in fchk.readlines()[line1:line2]:
			columns = line.split()
			for k in range(0, len(columns)):
				values.append(int(columns[k]))
		return values

	def fchk_get_first_double_values(self, phrase):
		values = []
		info = first_occurrence(self.fchkfile, phrase, 1)
		line1, nvalues = info[0], info[1]
		nlines = int(nvalues/5)
		if nlines * 5 < nvalues : nlines += 1
		line2 = line1 + nlines
		fchk = open(self.fchkfile, "r")
		#print 'phrase=', phrase, 'line1=', line1, 'nlines=', nlines
		for line in fchk.readlines()[line1:line2]:
			ncols = int((len(line)+1)/16)
			columns = []
			for k in range(0, ncols):
				columns.append(line[k*16:(k+1)*16])
			#columns = line.split()
			for k in range(0, len(columns)):
				values.append(float(columns[k]))
		return values

	def fchk_get_last_double_values(self, phrase):
		values = []
		info = last_occurrence(self.fchkfile, phrase, 1)
		line1, nvalues = info[0], info[1]
		nlines = int(nvalues/5)
		if nlines * 5 < nvalues : nlines += 1
		line2 = line1+nlines
		fchk = open(self.fchkfile, "r")
		#print('phrase=', phrase, 'line1=', line1, 'nlines=', nlines)
		for line in fchk.readlines()[line1:line2]:
			ncols = int((len(line)+1)/16)
			columns = []
			for k in range(0, ncols):
				columns.append(line[k*16:(k+1)*16])
			#columns = line.split()
			for k in range(0, len(columns)):
				values.append(float(columns[k]))
		return values

	def get_geometry(self):
		coords = self.fchk_get_last_double_values("Current cartesian coordinates")
		for k in range(0, 3*self.natoms): coords[k] *= BOHR
		atmnum = self.fchk_get_last_integer_values("Atomic numbers")
		self.geometry = Molecule()
		self.geometry.setup2(atmnum, coords)
		#self.geometry.mprint("read-in fchk geometry")

	def get_basis(self):
		self.nshl      =int(get_first_value(self.fchkfile, "Number of contracted shells", -1))
		self.nshl_prim =int(get_first_value(self.fchkfile, "Number of primitive shells", -1)) 
		types   = self.fchk_get_last_integer_values("Shell types")
		nprims  = self.fchk_get_last_integer_values("Number of primitives per shell")
		map     = self.fchk_get_last_integer_values("Shell to atom map")
		expo    = self.fchk_get_last_double_values("Primitive exponents")
		coef    = self.fchk_get_first_double_values("Contraction coefficients")
		coef_sp = self.fchk_get_last_double_values("P(S=P) Contraction coefficients")
		self.basis = BasisSet()
		self.basis.setup(self.geometry, types, nprims, map, expo, coef, coef_sp)
		"""
		print 'nbas=', self.basis.nbas
		print 'nbas_p=', self.basis.nbas_p
		print 'nshl=', self.basis.nshl
		print 'nshl_p=', self.basis.nshl_p
		for k in range(0, self.basis.nshl):
			print 'k=', k, 'coords=', self.basis.shells[k].coords
			print 'k=', k, 'expo=', self.basis.shells[k].expo
			print 'k=', k, 'coef=', self.basis.shells[k].coef
			if self.basis.shells[k].type == -1:
				print 'k=', k, 'coef_sp=', self.basis.shells[k].coef_sp
		"""

	def get_mos(self, type='mo'):
		if type == 'mo':
			return self.fchk_get_last_double_values("Alpha MO coefficients")
		elif type == 'b-mo':
			return self.fchk_get_last_double_values("Beta MO coefficients")
		elif type == 'frz':
			return self.fchk_get_last_double_values("Alpha FRZ coefficients")
		elif type == 'almo':
			return self.fchk_get_last_double_values("Alpha ALMO coefficients")
		elif type == 'rs':
			return self.fchk_get_last_double_values("Alpha RS coefficients")
		elif type == 'lowdin':
			return self.fchk_get_last_double_values("Alpha LOWDIN coefficients")
		elif type == 'proj':
			return self.fchk_get_last_double_values("Alpha Proj coefficients")
		elif type == 'pcanon':
			return self.fchk_get_last_double_values("Alpha pseudocanonical coeff")
		elif type == 'ibo-a':
			return self.fchk_get_last_double_values("Localized Alpha MO Coefficients (IBO)")
		elif type == 'ibo-b':
			return self.fchk_get_last_double_values("Localized Beta  MO Coefficients (IBO)")
		elif type == 'oslo-a':
			return self.fchk_get_last_double_values("Localized Alpha MO Coefficients (OSLO)")
		elif type == 'oslo-b':
			return self.fchk_get_last_double_values("Localized Beta  MO Coefficients (OSLO)")

	def get_mo_energies(self):
		return self.fchk_get_last_double_values("Alpha Orbital Energies")
	
	def get_nbos(self):
		return self.fchk_get_last_double_values("Alpha NBOs")
	
	def get_alpha_amplitudes(self):
		return self.fchk_get_last_double_values("Alpha Amplitudes")

	def get_square_matrix_from_uppertriangle(self, phrase):
		nbas = self.nbas
		array_upt = self.fchk_get_last_double_values(phrase)
		array = np.zeros(nbas*nbas)
		offset = 0 
		for j in range(0, nbas):
			for i in range(0, j+1):
				array[i+j*nbas] = array[j+i*nbas] = array_upt[offset+i]
			offset += j+1 
		return array
	    
	def get_overlap_matrix(self):
	        return self.get_square_matrix_from_uppertriangle("Overlap Matrix")
	
	def get_scf_density_matrix(self):
	        return self.get_square_matrix_from_uppertriangle("Total SCF Density")

	""" this two functions will be removed """
	def get_transition_density_matrix(self, istate):
		mos = self.get_mos()
		Co  = np.zeros(self.nbas*self.noa)
		Cv  = np.zeros(self.nbas*self.nva)
		for k in range(0, self.nbas*self.noa):
			Co[k] = mos[k]
		for k in range(0, self.nbas*self.nva):
			Cv[k] = mos[k+self.nbas*self.noa]
		X_all = self.get_alpha_amplitudes()
		X     = np.zeros(self.nva * self.noa)
		for k in range(0, self.nva*self.noa):
			X[k] = X_all[k+(istate-1)*self.nva*self.noa]
		return compute_transition_density_matrix(self.nbas, self.noa, self.nva, Co, Cv, X)
		
	def get_difference_density_matrix(self, ddm, attdm, detdm, istate):
		mos = self.get_mos()
		Co  = np.zeros(self.nbas*self.noa)
		Cv  = np.zeros(self.nbas*self.nva)
		for k in range(0, self.nbas*self.noa):
			Co[k] = mos[k]
		for k in range(0, self.nbas*self.nva):
			Cv[k] = mos[k+self.nbas*self.noa]
		X_all = self.get_alpha_amplitudes()
		X     = np.zeros(self.nva * self.noa)
		for k in range(0, self.nva*self.noa):
			X[k] = X_all[k+(istate-1)*self.nva*self.noa]
		compute_difference_density_matrix(ddm, attdm, detdm, self.nbas, self.noa, self.nva, Co, Cv, X)
		#matrix_print_2(ddm, self.nbas, self.nbas, 6, "ddm in fchk")
		return

	def reorder_mo_coef(self, C):
		norb = C.shape[0]
		print("norb:", norb)
		print("basis_order", self.basis.basis_order)
		C2 = np.zeros((norb, self.nbas))
		for i in range(0, norb):
			for j in range(0, self.nbas):
	        		C2[i,j] = C[i,self.basis.basis_order[j]]
		return C2

		
		

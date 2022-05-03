# -*- coding: utf-8 -*-
#
# Copyright (c), the QMAnalysis development team
#
#
#

""" basis set information"""

import numpy, sys, os
from .molecule import Molecule
from .util import *

def nMu(L):
	""" number of basis functions in a given angular-momentum shell """

	if L == 0:    return 1
	elif L == 1:  return 3
	elif L == -1: return 4
	elif L == 2:  return 6
	elif L == -2: return 5
	elif L == 3:  return 10
	elif L == -3: return 7

def nMu2(L):
	""" for separating grid values of basis functions into x, y, and z components """

	if L == 0:              return 1   
	elif L == 1 or L == -1: return 2 
	elif L == 2 or L == -2: return 3 
	elif L == 3 or L == -3: return 4

#basis function order is different in Q-Chem and FChk 
def order_for_a_shell(L):
	if L == 0:	return [0]
	elif L == 1: 	return [0,1,2]
	elif L == -1: 	return [0,1,2,3]
	elif L == 2: 	return [0,2,5,1,3,4]
	elif L == -2: 	return [2,3,1,4,0]
	elif L == 3: 	return [0,3,9,2,1,4,7,8,6,5]
	elif L == -3: 	return [3,4,2,5,1,6,0]

def basis_scale_factor(exponent, type, sp):
	""" normalization factors """

	Norm = pow(2.0/PI, 0.75);
	if type == 0 or (type == -1 and sp == 0):   Norm *= pow(exponent, 0.75)
	elif type == 1 or (type == -1 and sp == 1): Norm *= 2.0*pow(exponent, 1.25)
	elif type == 2 or type == -2:               Norm *= 4.0*pow(exponent, 1.75)
	elif type == 3 or type == -3:               Norm *= 8.0*pow(exponent, 2.25)
	
	return Norm

class Shell():
	""" a s, sp, p, d shell """
	""" f and higher functions not supported yet """

	def __init__(self):
		if 0 == 1: write("initialization a shell\n")

	def setup(self, coords, type, nprim, map, expo, coef, coef_sp):
		self.coords  = coords
		self.type    = type
		self.nprim   = nprim
		self.map     = map
		self.expo    = expo
		self.coef    = coef
		for k in range(0, nprim):
			self.coef[k] *= basis_scale_factor(expo[k], type, 0)
		if len(coef_sp) > 0:
			self.coef_sp = coef_sp
			if type == -1:
				for k in range(0, nprim):
					self.coef_sp[k] *= basis_scale_factor(expo[k], type, 1)

class BasisSet():
	""" Basis set information for a molecule """
	def __init__(self):
		if 0 == 1: write("initialization a shell\n")

	def setup(self, molecule, types, nprims, map, expo, coef, coef_sp):
		""" setup a basis for a molecule """

		self.nshl     = len(types)
		self.nshl_p   = len(expo)
		self.shells   = []
		self.offset   = []
		self.offset_p = []
		self.coords   = molecule.coords
		count = 0 
		count_basis   = 0
		count_basis_p = 0
		self.basis_order = []
		for k in range(0, self.nshl):
			nprim  = nprims[k]
			count2 = count + nprim
			shell = Shell()
			shell.setup(molecule.coords[3*map[k]-3:3*map[k]], types[k], nprims[k], map[k], expo[count:count2], coef[count:count2], coef_sp[count:count2])
			self.shells.append(shell)
			self.offset.append(count_basis)
			self.offset_p.append(count_basis_p)
			count += nprim
			order_shell = order_for_a_shell(types[k])
			for m in range(0, nMu(types[k])):
				self.basis_order.append(count_basis+order_shell[m])
			count_basis   += nMu(types[k])
			count_basis_p += nMu2(types[k])*nprim
		self.nbas   = count_basis
		self.nbas_p = count_basis_p


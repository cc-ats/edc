# -*- coding: utf-8 -*-
#
# Copyright (c), the QMAnalysis development team
#
#
#

""" compute basis function value on a user-defined grid """

import numpy, os, sys, time
import numpy as np
from numpy import linalg as LA
from sklearn.preprocessing import normalize
from math import exp
write = sys.stdout.write
from .molecule import Molecule
from .math_util import *
from .util import *

def project_onto_IC(vcc, B_matrix, Ginv):

	Bv = np.dot(vcc, B_matrix)
	#print("Bv=", Bv)
	vic= np.dot(Bv, Ginv)
	#print("vic=", vic)
	return vic
	
def get_G_inverse(B_matrix):
	nic = B_matrix.shape[1]
	print("nic:", nic)

	#G = B B^t
	G = np.dot(B_matrix.transpose(), B_matrix)
	matrix_print_2d(G, 6, "G")

	eval, evec = linalg.eigh(G)
	matrix_print_1d(eval, 1, nic, 6, "eval")
	matrix_print_2d(evec, 6, "evec")

	alpha = -1
	evec2 = np.zeros((nic,nic))
	for j in range(0, nic):
		if eval[j] < 0.00000001:
			scale = 0.0
		else: 
			scale = pow(eval[j], 0.5*alpha)
		for i in range(0, nic):	
			evec2[i,j] = evec[i,j] * scale
	Ginv = np.dot(evec2, evec2.transpose())
	matrix_print_2d(Ginv, 6, "Ginv")

	#make sure Ginv is the inverse of G
	prod = np.dot(G, Ginv)
	matrix_print_2d(prod, 6, "G * Ginv")
	
	return Ginv

def get_internal_coordinates(mol):

	natoms = mol.natoms

	nic = 0
	types = []
	ic_atoms = []
	B_matrix = np.zeros(3*natoms * (100*natoms))
	
	#bonds first
	for i in range(0, natoms):
		for j in range(i+1, natoms):
			if mol.connectivity[i+j*natoms] == 1:
				types.append(1)  
				ic_atoms.append([i,j])
				nic += 1

	#angles
	for i in range(0, natoms):
		for j in range(0, natoms):
			for k in range(j+1, natoms):
				if mol.connectivity[i+j*natoms] == 1 and mol.connectivity[i+k*natoms] == 1:
					types.append(2)
					ic_atoms.append([j, i, k]) 
					nic += 1 
	#dihedral angles
	#k - i - j - l
	for i in range(0, natoms):
		for j in range(i+1, natoms):
			if mol.connectivity[i+j*natoms] == 1:
				for k in range(0, natoms):
					if k != j and mol.connectivity[i+k*natoms] == 1:
						for l in range(0, natoms):
							if l != i and l != k and mol.connectivity[j+l*natoms] == 1:
								types.append(3)
								ic_atoms.append([k, i, j, l]) 
								nic += 1

	#B matrix
	B_matrix = np.zeros((3*natoms,nic))
	for ic in range(0, nic):
		type = types[ic]
		v = np.zeros(3*natoms)
		if type == 1:  #bond
			i = ic_atoms[ic][0]
			j = ic_atoms[ic][1]
			vij = np.subtract(mol.coords[3*i:3*i+3], mol.coords[3*j:3*j+3])
			for m in range(0, 3):
				vij[m] /= mol.distances[i+j*natoms]/BOHR #normalize
				v[3*i+m] =  vij[m]
				v[3*j+m] = -vij[m]
		elif type == 2:  #angle
			j = ic_atoms[ic][0]
			i = ic_atoms[ic][1]
			k = ic_atoms[ic][2]
			vij = np.subtract(mol.coords[3*i:3*i+3], mol.coords[3*j:3*j+3])
			vik = np.subtract(mol.coords[3*i:3*i+3], mol.coords[3*k:3*k+3])
			dij = mol.distances[i+j*natoms]/BOHR
			dik = mol.distances[i+k*natoms]/BOHR
			fj = np.zeros(3)
			fk = np.zeros(3)
			cos = 0.0
			for m in range(0, 3):
				cos += vij[m] * vik[m] / (dij * dik)
			if abs(cos) > 0.01:
				#I think the following is off by -1, but just to be consistent with Jon Baker's code
				s = 1.0 / sqrt(1 - cos*cos)
				for m in range(0, 3):
					fj[m] = s * (- cos * vij [m] / (dij * dij) + vik[m] / (dij * dik) )  
					fk[m] = s * (- cos * vik [m] / (dik * dik) + vij[m] / (dij * dik) )  
					v[3*j+m] = fj[m]
					v[3*k+m] = fk[m]
					v[3*i+m] -= fj[m] + fk[m]
			else: 
				print("need to fix this")
		elif type == 3:  #angle
			k = ic_atoms[ic][0]
			i = ic_atoms[ic][1]
			j = ic_atoms[ic][2]
			l = ic_atoms[ic][3]
			vik = np.subtract(mol.coords[3*i:3*i+3], mol.coords[3*k:3*k+3])
			vij = np.subtract(mol.coords[3*i:3*i+3], mol.coords[3*j:3*j+3])
			vjl = np.subtract(mol.coords[3*j:3*j+3], mol.coords[3*l:3*l+3])
			dik = mol.distances[i+k*natoms]/BOHR
			dij = mol.distances[i+j*natoms]/BOHR
			djl = mol.distances[j+l*natoms]/BOHR
			dotkij = np.dot(vik, vij)
			dotijl = np.dot(vij, vjl)
			crosskij = np.cross(vik, vij)
			crossijl = np.cross(vij, vjl)
			norm_kij = LA.norm(crosskij)
			norm_ijl = LA.norm(crossijl)
			cos = 0.0
			for m in range(0, 3):
				cos += vik[m] * vjl[m] 
			cos /= norm_kij * norm_ijl  #dihedral angle
			fik = np.zeros(3)
			fjl = np.zeros(3)
			if abs(cos) > 0.01:
				s = 1.0 / sqrt(1 - cos*cos)
				for m in range(0, 3):
					fik[m] =    s * dij * crosskij[m] / (norm_kij * norm_kij)
					fjl[m] =    s * dij * crossijl[m] / (norm_ijl * norm_ijl)
					v[3*k+m] =  fik[m]
					v[3*i+m] = -fik[m]
					v[3*j+m] = -fjl[m]
					v[3*l+m] =  fjl[m]
					v[3*i+m] += (fik[m] * dotkij + fjl[m] * dotijl) / (dij*dij)
					v[3*j+m] -= (fik[m] * dotkij + fjl[m] * dotijl) / (dij*dij)
			else: 
				print("need to fix this")

		for m in range(0, 3*natoms):
			B_matrix[m][ic] = v[m]

	print_internal_coordinates(nic, types, ic_atoms, mol)
	matrix_print_2d(B_matrix, 6, "B_matrix")

	return types, ic_atoms, B_matrix

def print_internal_coordinates(nic, types, ic_atoms, mol):
	type_names = ['bond', 'angle', 'dihedral']
	write("Internal Coordinates:\n")
	for ic in range(0, nic):
		write("%3d: %8s" % (ic+1, type_names[types[ic]-1]))
		for iat in range(0, len(ic_atoms[ic])):
			i = ic_atoms[ic][iat]
			if iat > 0: write(" - ")
			write("%3s%-3d" % (mol.atmsym[i],i+1))
		write("\n")
	write("\n")


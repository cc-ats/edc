# -*- coding: utf-8 -*-
#
# Copyright (c), the QMAnalysis development team
#
#
#

""" Parser for Q-Chem output file """

import numpy, sys, os
from .molecule import *

def check_connectivity(molecule, nsub):

	natoms = molecule.natoms
	distances = np.zeros(natoms*natoms)
	connectivity = []
	for k in range(0, natoms*natoms): connectivity.append(0)
	get_distance_and_connectivity(distances, connectivity, molecule)

	#make sure the last nsub atoms is connecting to the rest by one covalent bond
	nbonds = 0
	atom1  = -1
	atom2  = -1
	for i in range(0, natoms-nsub):
		for j in range(natoms-nsub, natoms):
			if connectivity[i+j*natoms] == 1:
				nbonds += 1
				atom1  = i 
				atom2  = j

	if nbonds == 0 or nbonds > 1:
		print(nbonds, " bonds found between main molecule and substituent, fix it")
		print_distance_and_connectivity(distances, connectivity, molecule)
		sys.exit()

	order = []
	for k in range(0, natoms):
		order.append(k)
	if atom1 != natoms-nsub-1:
		for k in range(atom1,natoms-nsub-1):
			order[k] = k+1
		order[natoms-nsub-1] = atom1
	if atom2 != natoms-nsub:
		order[natoms-nsub] = atom2
		for k in range(natoms-nsub+1, atom2+1):
			order[k] = k-1
	print("order=", order)
	
	xyz2 = np.zeros(3*natoms)
	atmnum2 = []
	for k in range(0, natoms):
		for m in range(0, 3):
			xyz2[m+3*k] = molecule.coords[m+3*order[k]] * BOHR
		atmnum2.append(molecule.atmnum[order[k]])

	molecule2 = Molecule()
	molecule2.setup2(atmnum2, xyz2)
	molecule2.set_charge_multiplicity(molecule.charge, molecule.multiplicity)
	molecule2.nfrag = 2
	molecule2.frag_charges = [molecule.charge+1,-1]
	molecule2.frag_index = []
	for k in range(0, natoms):
		if k < natoms-nsub:
			molecule2.frag_index.append(0)
		else:
			molecule2.frag_index.append(1)
	print("frag_index:", molecule2.frag_index)

	return molecule2

#replace second fragment with h
def replace_sub_with_h(molecule):

	atom1 = 0
	xyz2 = []
	atmnum2 = []
	for k in range(0, molecule.natoms):
		if molecule.frag_index[k] == 0:
			for m in range(0, 3): 
				xyz2.append(molecule.coords[m+3*k] * BOHR)
			atmnum2.append(molecule.atmnum[k])
			atom1 = k
	atom2 = atom1 + 1

	d2 = 0.0
	v = []
	for m in range(0, 3):
		v.append(molecule.coords[m+3*atom2] - molecule.coords[m+3*atom1])
		d2 += v[m]*v[m]
	distance = sqrt(d2) * BOHR
	scale = 1.086 / distance
	for m in range(0, 3):
		xyz2.append((molecule.coords[m+3*atom1] + scale*v[m])*BOHR)
	atmnum2.append(1)
	
	molecule2 = Molecule()
	molecule2.setup2(atmnum2, xyz2)
	molecule2.set_charge_multiplicity(molecule.charge, molecule.multiplicity)
	molecule2.nfrag = 2
	molecule2.frag_charges = [molecule.charge+1,-1]
	molecule2.frag_index = []
	for k in range(0,atom1+1):
		molecule2.frag_index.append(0)
	molecule2.frag_index.append(1)
	print("frag_index:", molecule2.frag_index)

	return molecule2

#replace first fragment with benzene
def replace_main_with_benzene(molecule):

	natoms = molecule.natoms
	distances = np.zeros(natoms*natoms)
	connectivity = []
	for k in range(0, natoms*natoms): connectivity.append(0)
	get_distance_and_connectivity(distances, connectivity, molecule)

	atom1 = 0
	for k in range(0, natoms):
		if molecule.frag_index[k] == 0:
			atom1 = k
	atom2 = atom1+1

	atom3 = 0	
	for k in range(0, atom1):
		if connectivity[k+atom1*natoms] == 1:
			atom3 = k
	v1 = []
	v2 = []
	for m in range(0, 3):
		v1.append(molecule.coords[m+3*atom1] - molecule.coords[m+3*atom2])
		v2.append(molecule.coords[m+3*atom1] - molecule.coords[m+3*atom3])
	d2 = np.dot(v1, v1)
	for m in range(0, 3):
		v1[m] *= 1/sqrt(d2)    # unit long
	d12 = np.dot(v1, v2)
	for m in range(0, 3):
		v2[m] -= d12 * v1[m]
	d2 = np.dot(v2, v2)
	for m in range(0, 3):
		v2[m] *= 1/sqrt(d2)    # unit long

	RCC = 1.400
	RCH = 1.086

	COS = sqrt(3.)/2

	natoms2 = natoms - (atom1 - 10)
	xyz2 = np.zeros(natoms2*3)
	atmnum2 = [6,1,6,1,6,1,6,1,6,1,6]
	for m in range(0, 3):
		xyz2[m+3*0]  = molecule.coords[m+3*atom1]*BOHR + RCC * 0.5 * v1[m] + RCC * COS * v2[m]
		xyz2[m+3*1]  = xyz2[m+3*0] - RCH * 0.5 * v1[m] + RCH * COS * v2[m]
		xyz2[m+3*2]  = molecule.coords[m+3*atom1]*BOHR + RCC * 0.5 * v1[m] - RCC * COS * v2[m]
		xyz2[m+3*3]  = xyz2[m+3*2] - RCH * 0.5 * v1[m] - RCH * COS * v2[m]
		xyz2[m+3*4]  = xyz2[m+3*0] + RCC * v1[m]
		xyz2[m+3*5]  = xyz2[m+3*4] + RCH * 0.5 * v1[m] + RCH * COS * v2[m]
		xyz2[m+3*6]  = xyz2[m+3*2] + RCC * v1[m]
		xyz2[m+3*7]  = xyz2[m+3*6] + RCH * 0.5 * v1[m] - RCH * COS * v2[m]
		xyz2[m+3*8]  = molecule.coords[m+3*atom1]*BOHR + RCC * 2.0 * v1[m]
		xyz2[m+3*9]  = xyz2[m+3*8] + RCH * v1[m] 
		xyz2[m+3*10] = molecule.coords[m+3*atom1]*BOHR

	for k in range(atom1+1, natoms):
		for m in range(0, 3):
			xyz2[m+3*(11+k-atom1-1)] = molecule.coords[m+3*k]*BOHR
		atmnum2.append(molecule.atmnum[k])
		
	molecule2 = Molecule()
	molecule2.setup2(atmnum2, xyz2)
	molecule2.nfrag = 2
	molecule2.frag_charges = [1,-1]
	molecule2.frag_index = []
	for k in range(0, 11):
		molecule2.frag_index.append(0)
	for k in range(atom1+1, natoms):
		molecule2.frag_index.append(1)
	print("frag_index:", molecule2.frag_index)
 		
	return molecule2




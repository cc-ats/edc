# -*- coding: utf-8 -*-
#
# Copyright (c), the QMAnalysis development team
#
#
#

""" setup cubic grid and compute mo and density value on the grid """


import numpy, os, sys, time
import numpy as np
from .molecule import Molecule
from .util import *
from scipy import sparse
from scipy import stats
write = sys.stdout.write

def mo_grid_data(mo, basis_on_grid, use_csr):

	""" compute grid data for a molecule orbital """
	time1 = time.time()
	nbasis = len(mo)
	ngrid  = int(len(basis_on_grid)/nbasis)
	#print('nbasis=', nbasis, 'ngrid=', ngrid)
	print("use_csr: ", use_csr)

	if use_csr != 1:
		basis = np.reshape(basis_on_grid, (-1, ngrid))
		grid_data = np.dot(mo, basis)
	else:
		basis = np.reshape(basis_on_grid, (-1, ngrid))
		basis_sp = sparse.csr_matrix(basis)
		mo_sp = sparse.csr_matrix(mo)
		grid_data_sp = mo_sp * basis_sp
		grid_data = np.ravel(grid_data_sp.toarray())
		print("basis_sp", basis_sp.shape, "mo_sp", mo_sp.shape, "grid_data_sp", grid_data_sp.shape)
		print("nonzeros", basis_sp.count_nonzero(), mo_sp.count_nonzero(), grid_data_sp.count_nonzero())

	time2 = time.time()
	write("time for compute mo on the grid:  %7.1f\n" % (time2 - time1))
	return grid_data

def den_grid_data(sm, basis_on_grid, use_csr):

	""" compute grid data associated with a square matrix """
	time1 = time.time()
	nbasis = len(sm[0])
	ngrid  = int(len(basis_on_grid)/nbasis)
	#print('nbasis=', nbasis, 'ngrid=', ngrid)

	if use_csr != 1:
		basis = (np.reshape(basis_on_grid, (-1, ngrid))).transpose()
		sm_basis = np.dot(basis, sm)
	else: 
		basis = (np.reshape(basis_on_grid, (-1, ngrid))).transpose()
		basis_sp = sparse.csr_matrix(basis)
		sm2 = np.zeros((nbasis,nbasis))
		sparsity_sm = 0
		for k in range(0, nbasis):
			for l in range(0, nbasis):
				if abs(sm[k,l]) > 0.000001:
					sm2[k,l] = sm[k,l]
					sparsity_sm += 1
				else:
					sm2[k,l] = 0.0
		sparsity_sm /= nbasis*nbasis
		print("sparsity_sm:", sparsity_sm)
		sm_sp = sparse.csr_matrix(sm2)
		sm_basis_sp = basis_sp * sm_sp
		sm_basis = sm_basis_sp.toarray()
		time2 = time.time()
		write("time for matrix multiply:  %7.1f\n" % (time2 - time1))
	grid_data = np.zeros(ngrid)
	for k in range(0, ngrid):
		grid_data[k] = np.dot(basis[k], sm_basis[k]) 

	time3 = time.time()
	write("time for compute density on the grid:  %7.1f\n" % (time3 - time1))
	return grid_data

def average_local_ionization_energy_grid_data(density_matrix, mos, orbital_energies, basis_on_grid, nbasis, nalpha, use_csr):
	ngrid = int(len(basis_on_grid)/nbasis)
	print('nbasis=', nbasis, 'nalpha=', nalpha, 'ngrid=', ngrid)

	#total density
	dm_total = np.reshape(density_matrix, (-1, nbasis))
	density_total = den_grid_data(dm_total, basis_on_grid, use_csr)

	#grid data
	grid_data = np.zeros(ngrid)
	for k in range(0, nalpha):
		mo_k = np.reshape(mos[k*nbasis:(k+1)*nbasis], (-1, nbasis))
		dm_k = np.dot(mo_k.transpose(), mo_k)
		print("k:", k, "dm shape", dm_k.shape, "orbital energy", abs(orbital_energies[k]))
		density_k = den_grid_data(dm_k, basis_on_grid, use_csr)
		for m in range(0, ngrid):
			grid_data[m] += 2 * abs(orbital_energies[k]) * 27.211 * density_k[m] / density_total[m]

	vmin = np.amin(grid_data)
	vmax = np.amax(grid_data)
	vave = np.average(grid_data)
	print("vmin=", vmin, "vmax=", vmax, "vave=", vave)
	#grid_data -= vmin
	
	return grid_data

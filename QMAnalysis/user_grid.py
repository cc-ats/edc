# -*- coding: utf-8 -*-
#
# Copyright (c), the QMAnalysis development team
#
#
#

""" compute basis function value on a user-defined grid """

import numpy, os, sys, time
import numpy as np
from math import exp
write = sys.stdout.write
from .molecule import Molecule
from .util import *
from .basis_set import *

def get_user_grid_data_for_basis(grid, basis, natoms):

	""" compute grid value for basis functions """

	nshl   = basis.nshl
	nbas   = basis.nbas
	nbas_p = basis.nbas_p
	#print 'nshl=', nshl, 'nbas=', nbas, 'nbas_p=', nbas_p

	ngrid = int(len(grid)/3)
	grid_data = np.zeros(ngrid*nbas)

	#for k in range(0, nshl):
	#	print('k:', k, 'atom', basis.shells[k].map, 'bas1', basis.offset[k])
	print("ngrid: ", ngrid, "natoms: ", natoms)

	time1 = time.time()
	for igrid in range(0, ngrid):
		for iAtom in range(0, natoms):
			x = grid[3*igrid+0] - basis.coords[3*iAtom+0]
			y = grid[3*igrid+1] - basis.coords[3*iAtom+1]
			z = grid[3*igrid+2] - basis.coords[3*iAtom+2]
			r2 = x * x + y * y + z * z
			for k in range(0, nshl):    
				katom  = basis.shells[k].map - 1
				if iAtom == katom:
					type   = basis.shells[k].type
					bas1 = basis.offset[k]
					rad = 0.0
					rad_sp = 0.0
					for l in range(0, basis.shells[k].nprim):
						coeff = basis.shells[k].coef[l]
						rad += coeff * exp(-basis.shells[k].expo[l]*r2)
						if type == -1:
							coeff_sp = basis.shells[k].coef_sp[l]
							rad_sp += coeff_sp * exp(-basis.shells[k].expo[l]*r2)
					#print(igrid, k, iAtom, rad)
					if type == 0:	
						grid_data[bas1*ngrid+igrid] = rad;
					elif type == 1:
						grid_data[(bas1+0)*ngrid+igrid] = x*rad;
						grid_data[(bas1+1)*ngrid+igrid] = y*rad;
						grid_data[(bas1+2)*ngrid+igrid] = z*rad;
					elif type == -1:
						grid_data[(bas1+0)*ngrid+igrid] = rad;
						grid_data[(bas1+1)*ngrid+igrid] = x*rad_sp;
						grid_data[(bas1+2)*ngrid+igrid] = y*rad_sp;
						grid_data[(bas1+3)*ngrid+igrid] = z*rad_sp;
					elif type == 2:
						grid_data[(bas1+0)*ngrid+igrid] = IRT3*x*x*rad;
						grid_data[(bas1+1)*ngrid+igrid] = IRT3*y*y*rad;
						grid_data[(bas1+2)*ngrid+igrid] = IRT3*z*z*rad;
						grid_data[(bas1+3)*ngrid+igrid] = x*y*rad;
						grid_data[(bas1+4)*ngrid+igrid] = x*z*rad;
						grid_data[(bas1+5)*ngrid+igrid] = y*z*rad;
					elif type == -2:
						grid_data[(bas1+0)*ngrid+igrid] = 0.5*IRT3*(2*z*z-x*x-y*y)*rad;
						grid_data[(bas1+1)*ngrid+igrid] = x*z*rad;
						grid_data[(bas1+2)*ngrid+igrid] = y*z*rad;
						grid_data[(bas1+3)*ngrid+igrid] = 0.5*(x*x-y*y)*rad;
						grid_data[(bas1+4)*ngrid+igrid] = x*y*rad;
			
	#matrix_print_2(grid_data, grid.ngrid, nbasis, 10, "new grid_data")
	time2 = time.time()
	write("time for assembling basis: %7.1f\n" % (time2 - time2))

	return grid_data


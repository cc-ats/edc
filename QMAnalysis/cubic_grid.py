# -*- coding: utf-8 -*-
#
# Copyright (c), the QMAnalysis development team
#
#
#


""" setup cubic grid and compute mo and density value on the grid """

import numpy, os, sys, time
import numpy as np
from math import exp
write = sys.stdout.write
from .molecule import Molecule
from .util import *
from .basis_set import *

class CubicGrid():
	""" Basis set information for a molecule """
	def __init__(self):
		if 0 == 1: write("initialization a cubic grid\n")

	def setup(self, molecule, headroom0, spacing0):
		""" setup a basis for a molecule """

		headroom = headroom0 / BOHR
		spacing  = spacing0  / BOHR

		coord = np.reshape(molecule.coords, (-1, 3)) 
		cmin = np.amin(coord, 0) - headroom
		cmax = np.amax(coord, 0) + headroom
		ni = [0, 0, 0]
		for k in range(0, 3): 
			cmin[k] = int(cmin[k]/spacing)*spacing
			cmax[k] = int(cmax[k]/spacing)*spacing
			ni[k]   = int((cmax[k]-cmin[k]+0.001)/spacing) + 1 

		self.ngrid = ni[0]*ni[1]*ni[2]
		self.nx = ni[0]
		self.ny = ni[1]
		self.nz = ni[2]
		self.x = np.linspace(cmin[0], cmax[0], ni[0])
		self.y = np.linspace(cmin[1], cmax[1], ni[1])
		self.z = np.linspace(cmin[2], cmax[2], ni[2])
		self.xinc = self.yinc = self.zinc = spacing
		
		write('Grid information in Angstroms\n')
		write('   x: %3d min: %7.3f max: %7.3f incr: %7.3f \n' % (self.nx, self.x[0]*BOHR, self.x[self.nx-1]*BOHR, self.xinc*BOHR))
		write('   y: %3d min: %7.3f max: %7.3f incr: %7.3f \n' % (self.ny, self.y[0]*BOHR, self.y[self.ny-1]*BOHR, self.yinc*BOHR))
		write('   z: %3d min: %7.3f max: %7.3f incr: %7.3f \n' % (self.nz, self.z[0]*BOHR, self.z[self.nz-1]*BOHR, self.zinc*BOHR))

		self.pos = np.zeros(3*self.ngrid)
		for ix in range(0, ni[0]):
			for iy in range(0, ni[1]):
				for iz in range(0, ni[2]):
					ixyz = iz + iy*ni[2] + ix*ni[1]*ni[2]
					self.pos[3*ixyz+0] = self.x[ix]
					self.pos[3*ixyz+1] = self.y[iy]
					self.pos[3*ixyz+2] = self.z[iz]

	#set up cubic grid with ni, cmin, and cinc
	def setup2(self, ni, cmin, cinc):
		self.ngrid = ni[0]*ni[1]*ni[2]
		self.nx = ni[0]
		self.ny = ni[1]
		self.nz = ni[2]
		self.xinc = cinc[0]
		self.yinc = cinc[1]
		self.zinc = cinc[2]
		cmax = np.zeros(3)
		cmax[0] = cmin[0] + (ni[0]-1) * cinc[0]
		cmax[1] = cmin[1] + (ni[1]-1) * cinc[1]
		cmax[2] = cmin[2] + (ni[2]-1) * cinc[2]
		self.x = np.linspace(cmin[0], cmax[0], ni[0])
		self.y = np.linspace(cmin[1], cmax[1], ni[1])
		self.z = np.linspace(cmin[2], cmax[2], ni[2])
		
		write('Grid information in Angstroms\n')
		write('   x: %3d min: %7.3f max: %7.3f incr: %7.3f \n' % (self.nx, self.x[0]*BOHR, self.x[self.nx-1]*BOHR, self.xinc*BOHR))
		write('   y: %3d min: %7.3f max: %7.3f incr: %7.3f \n' % (self.ny, self.y[0]*BOHR, self.y[self.ny-1]*BOHR, self.yinc*BOHR))
		write('   z: %3d min: %7.3f max: %7.3f incr: %7.3f \n' % (self.nz, self.z[0]*BOHR, self.z[self.nz-1]*BOHR, self.zinc*BOHR))

		self.pos = np.zeros(3*self.ngrid)
		for ix in range(0, ni[0]):
			for iy in range(0, ni[1]):
				for iz in range(0, ni[2]):
					ixyz = iz + iy*ni[2] + ix*ni[1]*ni[2]
					self.pos[3*ixyz+0] = self.x[ix]
					self.pos[3*ixyz+1] = self.y[iy]
					self.pos[3*ixyz+2] = self.z[iz]

	def get_cubic_grid_data_for_basis(self, basis):
		""" compute grid value for basis functions """

		nshl   = basis.nshl
		nbas   = basis.nbas
		nbas_p = basis.nbas_p
		#print 'nshl=', nshl, 'nbas=', nbas, 'nbas_p=', nbas_p

		ngrid = self.ngrid
		nx    = self.nx
		ny    = self.ny
		nz    = self.nz
		#print 'ngrid=', ngrid, nx, ny, nz

		""" compute x, y, and z components """
		time1 = time.time()
		x_components = np.zeros(nbas_p*nx)
		y_components = np.zeros(nbas_p*ny)
		z_components = np.zeros(nbas_p*nz)
		for k in range(0, nshl):    
			iatom  = basis.shells[k].map
			type   = basis.shells[k].type
			basis1 = basis.offset_p[k]
			for i in range(0, nx):
				v = self.x[i] - basis.coords[3*iatom-3+0]
				for l in range(0, basis.shells[k].nprim):
					rad = exp(-basis.shells[k].expo[l]*v*v)
					x_components[(basis1+nMu2(type)*l)*nx+i] = rad 
					if abs(type) > 0:
						x_components[(basis1+1+nMu2(type)*l)*nx+i] = rad * v 
					if abs(type) > 1:
						x_components[(basis1+2+nMu2(type)*l)*nx+i] = rad * v * v 
			for i in range(0, ny):
				v = self.y[i] - basis.coords[3*iatom-3+1]
				for l in range(0, basis.shells[k].nprim):
					rad = exp(-basis.shells[k].expo[l]*v*v)
					y_components[(basis1+nMu2(type)*l)*ny+i] = rad 
					if abs(type) > 0: y_components[(basis1+1+nMu2(type)*l)*ny+i] = rad * v 
					if abs(type) > 1: y_components[(basis1+2+nMu2(type)*l)*ny+i] = rad * v * v 
			for i in range(0, nz):
				v = self.z[i] - basis.coords[3*iatom-3+2]
				for l in range(0, basis.shells[k].nprim):
					rad = exp(-basis.shells[k].expo[l]*v*v)
					z_components[(basis1+nMu2(type)*l)*nz+i] = rad 
					if abs(type) > 0: z_components[(basis1+1+nMu2(type)*l)*nz+i] = rad * v 
					if abs(type) > 1: z_components[(basis1+2+nMu2(type)*l)*nz+i] = rad * v * v 
	
		#matrix_print_2(x_components, nx, nbas_p, 6, "new x_component")
		#matrix_print_2(y_components, ny, nbas_p, 6, "new y_component")
		#matrix_print_2(z_components, nz, nbas_p, 6, "new z_component")
	
		time2 = time.time()
		write("time for computing x,y,z components: %7.1f\n" % (time2 - time1))
	
		grid_data = np.zeros(ngrid*nbas)
		for k in range(0, nshl):
			type = basis.shells[k].type
			bas1 = basis.offset[k]
			for l in range(0, basis.shells[k].nprim):
				coeff = basis.shells[k].coef[l]
				off   = basis.offset_p[k] + nMu2(type)*l
				off1, off2, off3 = off+1, off+2, off+3
				if type == 0:
					x0y0   = np.outer(x_components[off*nx:off1*nx], y_components[off*ny:off1*ny])
					x0y0z0 = np.outer(x0y0, z_components[off*nz:off1*nz])
					grid_data[bas1*ngrid:(bas1+1)*ngrid] += coeff*x0y0z0.ravel()
				elif type == 1:
					x0y0   = np.outer(x_components[off*nx:off1*nx],  y_components[off*ny:off1*ny])
					x1y0   = np.outer(x_components[off1*nx:off2*nx], y_components[off*ny:off1*ny])
					x0y1   = np.outer(x_components[off*nx:off1*nx],  y_components[off1*ny:off2*ny])
					x1y0z0 = np.outer(x1y0, z_components[off*nz:off1*nz])
					x0y1z0 = np.outer(x0y1, z_components[off*nz:off1*nz])
					x0y0z1 = np.outer(x0y0, z_components[off1*nz:off2*nz])
					grid_data[(bas1+0)*ngrid:(bas1+1)*ngrid] += coeff*x1y0z0.ravel()
					grid_data[(bas1+1)*ngrid:(bas1+2)*ngrid] += coeff*x0y1z0.ravel()
					grid_data[(bas1+2)*ngrid:(bas1+3)*ngrid] += coeff*x0y0z1.ravel()
				elif type == -1:
					coeff_sp = basis.shells[k].coef_sp[l]
					x0y0   = np.outer(x_components[off*nx:off1*nx],  y_components[off*ny:off1*ny])
					x1y0   = np.outer(x_components[off1*nx:off2*nx], y_components[off*ny:off1*ny])
					x0y1   = np.outer(x_components[off*nx:off1*nx],  y_components[off1*ny:off2*ny])
					x0y0z0 = np.outer(x0y0, z_components[off*nz:off1*nz])
					x1y0z0 = np.outer(x1y0, z_components[off*nz:off1*nz])
					x0y1z0 = np.outer(x0y1, z_components[off*nz:off1*nz])
					x0y0z1 = np.outer(x0y0, z_components[off1*nz:off2*nz])
					grid_data[(bas1+0)*ngrid:(bas1+1)*ngrid] += coeff*x0y0z0.ravel()
					grid_data[(bas1+1)*ngrid:(bas1+2)*ngrid] += coeff_sp*x1y0z0.ravel()
					grid_data[(bas1+2)*ngrid:(bas1+3)*ngrid] += coeff_sp*x0y1z0.ravel()
					grid_data[(bas1+3)*ngrid:(bas1+4)*ngrid] += coeff_sp*x0y0z1.ravel()
				elif type == 2:
					x0y0   = np.outer(x_components[off*nx:off1*nx],  y_components[off*ny:off1*ny])
					x1y0   = np.outer(x_components[off1*nx:off2*nx], y_components[off*ny:off1*ny])
					x2y0   = np.outer(x_components[off2*nx:off3*nx], y_components[off*ny:off1*ny])
					x0y1   = np.outer(x_components[off*nx:off1*nx],  y_components[off1*ny:off2*ny])
					x1y1   = np.outer(x_components[off1*nx:off2*nx], y_components[off1*ny:off2*ny])
					x0y2   = np.outer(x_components[off*nx:off1*nx],  y_components[off2*ny:off3*ny])
					x2y0z0 = np.outer(x2y0, z_components[off*nz:off1*nz])
					x1y1z0 = np.outer(x1y1, z_components[off*nz:off1*nz])
					x0y2z0 = np.outer(x0y2, z_components[off*nz:off1*nz])
					x1y0z1 = np.outer(x1y0, z_components[off1*nz:off2*nz])
					x0y1z1 = np.outer(x0y1, z_components[off1*nz:off2*nz])
					x0y0z2 = np.outer(x0y0, z_components[off2*nz:off3*nz])
					grid_data[(bas1+0)*ngrid:(bas1+1)*ngrid] += (coeff*IRT3)*x2y0z0.ravel()
					grid_data[(bas1+1)*ngrid:(bas1+2)*ngrid] += (coeff*IRT3)*x0y2z0.ravel()
					grid_data[(bas1+2)*ngrid:(bas1+3)*ngrid] += (coeff*IRT3)*x0y0z2.ravel()
					grid_data[(bas1+3)*ngrid:(bas1+4)*ngrid] += coeff*x1y1z0.ravel()
					grid_data[(bas1+4)*ngrid:(bas1+5)*ngrid] += coeff*x1y0z1.ravel()
					grid_data[(bas1+5)*ngrid:(bas1+6)*ngrid] += coeff*x0y1z1.ravel()
				elif type == -2:
					x0y0   = np.outer(x_components[off*nx:off1*nx],  y_components[off*ny:off1*ny])
					x1y0   = np.outer(x_components[off1*nx:off2*nx], y_components[off*ny:off1*ny])
					x2y0   = np.outer(x_components[off2*nx:off3*nx], y_components[off*ny:off1*ny])
					x0y1   = np.outer(x_components[off*nx:off1*nx],  y_components[off1*ny:off2*ny])
					x1y1   = np.outer(x_components[off1*nx:off2*nx], y_components[off1*ny:off2*ny])
					x0y2   = np.outer(x_components[off*nx:off1*nx],  y_components[off2*ny:off3*ny])
					x2y0z0 = np.outer(x2y0, z_components[off*nz:off1*nz])
					x1y1z0 = np.outer(x1y1, z_components[off*nz:off1*nz])
					x0y2z0 = np.outer(x0y2, z_components[off*nz:off1*nz])
					x1y0z1 = np.outer(x1y0, z_components[off1*nz:off2*nz])
					x0y1z1 = np.outer(x0y1, z_components[off1*nz:off2*nz])
					x0y0z2 = np.outer(x0y0, z_components[off2*nz:off3*nz])
					grid_data[(bas1+0)*ngrid:(bas1+1)*ngrid] += (coeff*0.5*IRT3)*(2*x0y0z2.ravel()-x2y0z0.ravel()-x0y2z0.ravel())
					grid_data[(bas1+1)*ngrid:(bas1+2)*ngrid] += coeff*x1y0z1.ravel()
					grid_data[(bas1+2)*ngrid:(bas1+3)*ngrid] += coeff*x0y1z1.ravel()
					grid_data[(bas1+3)*ngrid:(bas1+4)*ngrid] += (coeff*0.5)*(x2y0z0.ravel()-x0y2z0.ravel())
					grid_data[(bas1+4)*ngrid:(bas1+5)*ngrid] += coeff*x1y1z0.ravel()
	
		#matrix_print_2(grid_data, grid.ngrid, nbasis, 10, "new grid_data")
		time3 = time.time()
		write("time for assembling basis: %7.1f\n" % (time3 - time2))

		return grid_data



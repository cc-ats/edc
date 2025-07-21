# -*- coding: utf-8 -*-
#
# Copyright (c), the QMAnalysis development team
#
#
#

""" Parser for Q-Chem output file """
from numpy import ma
import numpy, sys, os
from .compute_grid_data import *
from .connectivity import *
from .cubic_grid import *
from .cube_file import *
from .make_plot import *
from .qchem_fchk_parser import *
from .user_grid import *
from .util import *
write = sys.stdout.write
from scipy import sparse
from scipy import stats
from scipy import linalg

class QMSystem():
	""" Q-Chem output file. """

	def __init__(self):
		#super(QChem, self).__init__(logname="QChem")
		#doing nothing
		if 1 == 0: write("initialization")

	def setup_from_qchem_fchk(self, fchkfile):
		"""setup molecule with atomic symbols"""
		self.fchkfile           = fchkfile
		self.fchk               = QChemFChk(fchkfile)
		self.system_name        = self.fchkfile[0:-4]
		if self.system_name[-4:-1] == '.in':
			self.system_name = self.system_name[0:-3]
		elif self.system_name[-5:-1] == '.inp':
			self.system_name = self.system_name[0:-4]

		#geometry
		self.molecule           = self.fchk.geometry
		self.charge             = self.molecule.charge
		self.multiplicity       = self.molecule.multiplicity
		self.natoms             = self.fchk.natoms
		self.nalpha             = self.fchk.nalpha
		self.nbasis             = self.fchk.nbas
		self.norb               = self.fchk.norb

		#basis, MOs, and density
		self.atomic_basis         = self.fchk.basis
		self.orbital_energies     = self.fchk.get_mo_energies()  
		self.molecular_orbitals   = self.fchk.get_mos()  
		self.total_density_matrix = self.fchk.get_scf_density_matrix()

		#overlap matrix
		#self.overlap_matrix       = self.fchk.get_overlap_matrix()

		#amplitudes from CIS/TDA calculations
		if first_occurrence(self.fchkfile, "Alpha Amplitudes") != -1:
			self.valpha       = self.norb - self.nalpha
			self.nova         = self.valpha*self.nalpha
			self.amplitudes   = self.fchk.get_alpha_amplitudes()
			self.nstates      = len(self.amplitudes)/self.nova
			write("nstates: %d\n" % self.nstates)

	def user_defined_grid(self, grid):
		ngrid = int(len(grid)/3)
		self.user_defined_grid = 1
		self.grid = np.zeros(ngrid*3)
		for k in range(0, ngrid*3): self.grid[k] = grid[k]

	""" spread one MO-like orbital over the grid """
	""" should work for both canonical orbitals and natural orbitals """

	def sparsity_analysis(self):
		def mythreshold(a, threshmin=None, threshmax=None, newval=0):
			a = ma.array(a , copy=True)
			mask = np.zeros(a.shape, dtype=bool)
			if threshmin is not None:
				mask |= (a < threshmin).filled(False)
			if threshmax is not None:
				mask |= (a > threshmax).filled(False)
			a[mask] = newval
			return a
		time1 = time.time()
		basis_on_grid2 = mythreshold(self.basis_on_grid, -0.000001, 0.000001)
		self.basis_on_grid -= basis_on_grid2
		ngrid = int(len(self.basis_on_grid)/self.nbasis)
		basis = np.reshape(self.basis_on_grid, (-1, ngrid))
		basis_sp = sparse.csr_matrix(basis)
		sparsity = 1.0 - 1.0 * basis_sp.count_nonzero() / ngrid / self.nbasis
		#sparsity = 0.0
		print("sparsity of basis_on_grid:", sparsity)
		self.use_csr = 0 
		if sparsity > 0.4 : self.use_csr = 1 
		time2 = time.time()
		write("time for sparsity analysis on basis_on_grid: %7.1f\n" % (time2 - time1))

	def plot_orbital_energies(self, plotfile, norbs=[10,10]):
		#print 10 occupied and 10 virtuals by default
		#print("norbs: ", norbs)
		if self.nalpha < norbs[0]:
			norbs[0] = self.nalpha
		if self.norb - self.nalpha < norbs[1]:
			norbs[1] = self.norb - self.nalpha

		norb = norbs[0] + norbs[1]
		energies = []
		occupancies = []
		for k in range(0, norb): 
			energies.append(self.orbital_energies[self.nalpha-norbs[0]+k])
			if (k < norbs[1]) :  occupancies.append(2)
			else: occupancies.append(0)
		print("energies: ", energies)
		print("occupancies: ", occupancies)
		make_orbital_energy_plot(plotfile, energies, occupancies)

	def plot_one_orbital(self, cube_file, img_file, mo_coefs, grid_size=0.1):
		#setup grid if necessary
		if not hasattr(self, 'grid'):
			self.grid = CubicGrid()		
			self.grid.setup(self.molecule, 2.0, grid_size)

		#generate basis data, if necessary
		if not hasattr(self, 'basis_on_grid'):
			#write("spread the basis over the grid\n")
			self.basis_on_grid = self.grid.get_cubic_grid_data_for_basis(self.atomic_basis)		
			self.use_scr = self.sparsity_analysis()

		#spread the MOs on the grid
		coefs = np.zeros(self.nbasis)
		scale = 0
		for k in range(0, self.nbasis):
			if scale == 0 and fabs(mo_coefs[k]) > 0.05:
				if mo_coefs[k] > 0:
					scale = 1
				else:
					scale = -1
		if scale == -1: print(" **** changing the sign of the orbital ****")
		for k in range(0, self.nbasis):
			coefs[k] = scale * mo_coefs[k]
		#print("coefs=", coefs)

		#mo_on_grid = mo_grid_data(mo_coefs, self.basis_on_grid, self.use_csr)
		mo_on_grid = mo_grid_data(coefs, self.basis_on_grid, self.use_csr)

		#generate cube file
		write_cube(cube_file, self.molecule, self.grid, mo_on_grid)

		#generate image file
		if img_file is not None:
			make_jmol_isosurface_plot(img_file, cube_file)


	def plot_selected_mos(self, mo_list, grid_size=0.1, type='mo'):
		if type == 'frz':
			self.frz_orbitals   = self.fchk.get_mos('frz')  
		elif type == 'almo':
			self.almo_orbitals   = self.fchk.get_mos('almo')  
		elif type == 'rs':
			self.rs_orbitals   = self.fchk.get_mos('rs')  
		elif type == 'lowdin':
			self.lowdin_orbitals   = self.fchk.get_mos('lowdin')  
		elif type == 'proj':
			self.proj_orbitals   = self.fchk.get_mos('proj')  
		elif type == 'pcanon':
			self.pcanon_orbitals   = self.fchk.get_mos('pcanon')  
		for k in range(0, len(mo_list)):

			#choose filenames for the cube and image data
			imo = mo_list[k]
			if imo <= self.nalpha:
				mo_name = self.system_name + "homo"
				if imo <= self.nalpha-1:
					mo_name += "-"+str(self.nalpha-imo)
			else: 
				mo_name = self.system_name + "lumo"
				if imo >= self.nalpha+2:
					mo_name += "+"+str(imo-self.nalpha-1)
			if type == 'frz' or type == 'almo' or type == 'rs' or type == 'lowdin' or type == 'proj' or type == 'pcanon':
				mo_name += "_" +type
			cube_file = mo_name+".cube"
			img_file  = mo_name+".png"

			#obtain mo coefficents
			if type == 'mo':
				mo_coefs  = self.molecular_orbitals[(imo-1)*self.nbasis:imo*self.nbasis]
			elif type == 'frz':
				mo_coefs  = self.frz_orbitals[(imo-1)*self.nbasis:imo*self.nbasis]
			elif type == 'almo':
				mo_coefs  = self.almo_orbitals[(imo-1)*self.nbasis:imo*self.nbasis]
			elif type == 'rs':
				mo_coefs  = self.rs_orbitals[(imo-1)*self.nbasis:imo*self.nbasis]
			elif type == 'lowdin':
				mo_coefs  = self.lowdin_orbitals[(imo-1)*self.nbasis:imo*self.nbasis]
			elif type == 'proj':
				mo_coefs  = self.proj_orbitals[(imo-1)*self.nbasis:imo*self.nbasis]
			elif type == 'pcanon':
				mo_coefs  = self.pcanon_orbitals[(imo-1)*self.nbasis:imo*self.nbasis]

			#generate cube and plot 
			self.plot_one_orbital(cube_file, img_file, mo_coefs, grid_size)

	""" spread density over the grid for a density-like matrix """
	""" should also work for a potential """
	def plot_one_density(self, cube_file, img_file, density_matrix, grid_size=0.1):
		#setup grid if necessary
		if not hasattr(self, 'grid'):
			self.grid = CubicGrid()		
			self.grid.setup(self.molecule, 2.0, grid_size)

		#generate basis data, if necessary
		if not hasattr(self, 'basis_on_grid'):
			#write("spread the basis over the grid\n")
			if self.user_defined_grid != 1:
				self.basis_on_grid = self.grid.get_cubic_grid_data_for_basis(self.atomic_basis)
				self.use_csr = self.sparsity_analysis()
			else:
				self.basis_on_grid = get_user_grid_data_for_basis(self.grid, self.atomic_basis, self.natoms)
				self.use_csr = 0

		#spread the total density over the grid
		#   using total density matrix
		dm = np.reshape(density_matrix, (-1, self.nbasis))
		self.density_on_grid = den_grid_data(dm, self.basis_on_grid, self.use_csr)
		if self.user_defined_grid == 1: return

		#generate cube files
		write_cube(cube_file, self.molecule, self.grid, self.density_on_grid)

		#find maximum value
		value = np.amax(np.absolute(self.density_on_grid))
		if value < 1.0: cutoff = 0.001
		else: cutoff = 0.02
		write("value: %f cutoff: %f\n" % (value, cutoff))

		#generate image file
		if img_file is not None:
			make_jmol_isosurface_plot(img_file, cube_file, cutoff)

	def plot_almo_density(self, type, name, ifrag, grid_size=0.1):
		if type == 'almo-frz':
			dims = self.fchk.fchk_get_last_integer_values("ALMO Data")
			print("dims: ", dims)
			if ifrag == 1:
				occ1, occ2 = 0, dims[2]
			else:
				occ1, occ2 = dims[2], dims[2]+dims[3]
			print("occ1=", occ1, "occ2=", occ2)
			frz_orbitals   = np.reshape(self.fchk.get_mos('frz')[self.nbasis*occ1:self.nbasis*occ2], (-1, self.nbasis))
			frz_density = 2.0 * np.dot(frz_orbitals.transpose(), frz_orbitals)
			almo_orbitals  = np.reshape(self.fchk.get_mos('almo')[self.nbasis*occ1:self.nbasis*occ2], (-1, self.nbasis))
			almo_density = 2.0 * np.dot(almo_orbitals.transpose(), almo_orbitals)
			diff_density = np.subtract(almo_density, frz_density).ravel()
			cube_file = name+"."+type+".frag"+str(ifrag)+".cube"
			img_file = name+"."+type+".frag"+str(ifrag)+".png"
			self.plot_one_density(cube_file, img_file,  diff_density, grid_size)
			
	""" compute the difference between qmmm and frz densities """
	def plot_qmmm_frz_density(self, type, name, ifrag, dm, grid_size=0.1):
		if type == 'qmmm-frz':
			dims = self.fchk.fchk_get_last_integer_values("ALMO Data")
			print("dims: ", dims)
			if ifrag == 1:
				occ1, occ2 = 0, dims[2]
			else:
				occ1, occ2 = dims[2], dims[2]+dims[3]
			print("occ1=", occ1, "occ2=", occ2)
			frz_orbitals   = np.reshape(self.fchk.get_mos('frz')[self.nbasis*occ1:self.nbasis*occ2], (-1, self.nbasis))
			frz_density = 2.0 * np.dot(frz_orbitals.transpose(), frz_orbitals)
			nbas1, nbas2, nbas = dims[0], dims[1], dims[0]+dims[1]
			print("nbas1=", nbas1, "nbas2=", nbas2, "nbas=", nbas)
			qmmm_density = np.zeros(nbas*nbas)
			for j in range(0, nbas2):
				for i in range(0, nbas2):
					qmmm_density[(nbas1+i)+(nbas1+j)*nbas] = dm[i+j*nbas2]
			diff_density = np.subtract(qmmm_density, frz_density.ravel())
			cube_file = name+"."+type+".frag"+str(ifrag)+".cube"
			img_file = name+"."+type+".frag"+str(ifrag)+".png"
			self.plot_one_density(cube_file, img_file,  diff_density, grid_size)

	""" spread the total density over the grid using total density matrix """
	def plot_total_density(self, grid_size=0.1):

		#choose filename for cube and image data
		den_name = self.system_name + 'total_density'
		cube_file = den_name+".cube"
		img_file  = den_name+".png"

		#make cube and image files
		self.plot_one_density(cube_file, img_file, self.total_density_matrix, grid_size)

	""" plot density (or mos) from input grid value """
	def plot_grid_data(self, filename, grid_data):

		cube_file = filename+".cube"
		img_file  = filename+".png"

		#generate cube files
		write_cube(cube_file, self.molecule, self.grid, grid_data)

		#find maximum value
		value = np.amax(np.absolute(grid_data))
		if value < 1.0: cutoff = 0.001
		else: cutoff = 0.02
		write("value: %f cutoff: %f\n" % (value, cutoff))

		#generate image file
		if img_file is not None:
			make_jmol_isosurface_plot(img_file, cube_file, cutoff)

	""" plot the average local ionization energy, Politzer et al, Can. J. Chem. 68, 1440 (1990)"""
	def plot_average_local_ionization_energy(self, grid_size=0.1):

		#choose filename for cube and image data
		den_name = self.system_name + 'average_local_ionization_energy'
		cube_file = den_name+".cube"
		img_file  = den_name+".png"

		#setup grid if necessary
		if not hasattr(self, 'grid'):
			self.grid = CubicGrid()		
			self.grid.setup(self.molecule, 2.0, grid_size)

		#generate basis data, if necessary
		if not hasattr(self, 'basis_on_grid'):
			#write("spread the basis over the grid\n")
			self.basis_on_grid = self.grid.get_cubic_grid_data_for_basis(self.atomic_basis)		
			self.use_scr = self.sparsity_analysis()

		#spread the total density over the grid
		value_on_grid = average_local_ionization_energy_grid_data(self.total_density_matrix, self.molecular_orbitals, self.orbital_energies, self.basis_on_grid, self.nbasis, self.nalpha, self.use_csr)

		#generate cube files
		write_cube(cube_file, self.molecule, self.grid, value_on_grid)

		#find maximum value
		value = np.amax(np.absolute(value_on_grid))
		if value < 1.0: cutoff = 0.001
		else: cutoff = 0.02
		write("value: %f cutoff: %f\n" % (value, cutoff))

		#generate image file
		if img_file is not None:
			make_jmol_isosurface_plot(img_file, cube_file, cutoff)

	def plot_tda_excited_state_densities(self, state_list, grid_size=0.1):

		for k in range(0, len(state_list)):

			#grab amplitudes and MO coefficients
			istate = state_list[k]
			X      = self.amplitudes[self.nova*(istate-1):self.nova*istate]
			Co     = self.molecular_orbitals[0:self.nbasis*self.nalpha]
			Cv     = self.molecular_orbitals[self.nbasis*self.nalpha:self.nbasis*self.norb]
			#print("Co length: ", len(Co), len(Cv), len(X))

			#compute transition, difference, attachment and detachment density matrices
			tdm    = compute_transition_density_matrix(self.nbasis, self.nalpha, self.valpha, Co, Cv, X)
			ddm, attdm, detdm  = [], [], []
			compute_difference_density_matrix(ddm, attdm, detdm, self.nbasis, self.nalpha, self.valpha, Co, Cv, X)

			#transition density
			den_name = self.system_name+'state'+str(istate)+'.transition_density'
			cube_file, img_file = den_name+".cube", den_name+".png"
			self.plot_one_density(cube_file, img_file, tdm, grid_size)
			
			#difference density
			den_name = self.system_name+'state'+str(istate)+'.difference_density'
			cube_file, img_file = den_name+".cube", den_name+".png"
			self.plot_one_density(cube_file, img_file, ddm)
			
			#attachment density
			den_name = self.system_name+'state'+str(istate)+'.attachment_density'
			cube_file, img_file = den_name+".cube", den_name+".png"
			self.plot_one_density(cube_file, img_file, attdm)
			
			#detachment density
			den_name = self.system_name+'state'+str(istate)+'.detachment_density'
			cube_file, img_file = den_name+".cube", den_name+".png"
			self.plot_one_density(cube_file, img_file, detdm)


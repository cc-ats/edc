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
from .math_util import *
from .molecule import Molecule


class QChemOut():
	""" Q-Chem output file. """

	def __init__(self, outputfile):
		self.outputfile = outputfile
		self.setup()

	def setup(self):
		self.has_optimized_geometry = -1
		self.has_last_geometry      = -1
		self.optimized_energy       = -100000.0
		opt = first_occurrence(self.outputfile, "Optimization Cycle:")
		if opt == -1:
			self.single_point()
		else:
			if first_occurrence(self.outputfile, "OPTIMIZATION CONVERGED") != -1:
				self.geometry_optimization()
			else:
				self.failed_geometry_optimization()

	def single_point(self):
		self.get_natoms()
		self.get_nbas()
		self.get_nele()
		self.get_input_geometry()
		self.get_scf_energy()
		#self.get_charges()

	def geometry_optimization(self):
		self.get_natoms()
		self.get_nbas()
		self.get_nele()
		self.get_input_geometry()
		self.get_optimized_geometry()
		self.get_scf_energy()
		self.get_optimized_energy()
		self.get_charges()

	def failed_geometry_optimization(self):
		self.get_natoms()
		self.get_nbas()
		self.get_nele()
		self.get_input_geometry()
		self.get_last_geometry()
		self.get_lowest_energy_geometry()
		self.get_scf_energy()
		self.get_charges()

	#There are 16 shells and 38 basis functions
	def get_nbas(self):
		self.nshl = int(get_first_value_2(self.outputfile, "basis functions", "shells", 2))
		self.nbas = int(get_first_value_2(self.outputfile, "basis functions", "shells", 5))

	#There are        8 alpha and        8 beta electrons
	def get_nele(self):
		self.nalpha = int(get_first_value_2(self.outputfile, "alpha", "beta", 2))
		self.nbeta  = int(get_first_value_2(self.outputfile, "alpha", "beta", 5))

	def get_natoms(self):
		line1 = first_occurrence(self.outputfile, "Standard Nuclear Orientation")
		line2 = first_occurrence(self.outputfile, "Nuclear Repulsion Energy")
		line3 = first_occurrence(self.outputfile, "Molecular Point Group")
		if line3 != -1:
			line2 = line3
		self.natoms = line2 - line1 - 4

	def get_input_geometry(self):
		line1 = first_occurrence(self.outputfile, "Standard Nuclear Orientation")
		outf = open(self.outputfile, "r", encoding = "ISO-8859-1")
		atmsym = []
		coords = []
		#sys.exit()
		for line in outf.readlines()[line1+2:line1+self.natoms+2]:
			columns = line.split()
			atmsym.append(columns[1])
			coords.append(float(columns[2]))
			coords.append(float(columns[3]))
			coords.append(float(columns[4]))
		self.input_geometry = Molecule()
		self.input_geometry.setup(atmsym, coords)

	def get_optimized_geometry(self):
		line1 = first_occurrence(self.outputfile, "OPTIMIZATION CONVERGED")
		outf = open(self.outputfile, "r",   encoding = "ISO-8859-1")
		atmsym = []
		coords = []
		for line in outf.readlines()[line1+4:line1+self.natoms+4]:
			columns = line.split()
			atmsym.append(columns[1])
			coords.append(float(columns[2]))
			coords.append(float(columns[3]))
			coords.append(float(columns[4]))
		if len(atmsym) == self.natoms and len(coords) == 3*self.natoms:
			self.optimized_geometry = Molecule()
			self.optimized_geometry.setup(atmsym, coords)
			self.has_optimized_geometry = 1

	def get_last_geometry(self):
		line1 = last_occurrence(self.outputfile, "Standard Nuclear Orientation")
		outf = open(self.outputfile, "r")
		atmsym = []
		coords = []
		for line in outf.readlines()[line1+2:line1+self.natoms+2]:
			columns = line.split()
			atmsym.append(columns[1])
			coords.append(float(columns[2]))
			coords.append(float(columns[3]))
			coords.append(float(columns[4]))
		if len(atmsym) == self.natoms and len(coords) == 3*self.natoms:
			self.last_geometry = Molecule()
			self.last_geometry.setup(atmsym, coords)
			self.has_last_geometry = 1

	def get_lowest_energy_geometry(self):
		index = []
		energies = get_all_values_with_index(self.outputfile, "Energy is", -1, index)
		lowest_energy = 0.0
		line1 = 0
		index_k = -1
		for k in range(0, len(energies)):
			if float(energies[k]) < lowest_energy:
				lowest_energy = float(energies[k])
				index_k = k+1
				line1 = index[k]
		line1 -= self.natoms + 4  # step back natoms+4 lines
		print('Lowest Energy: %f is found for the %d-th geometry' % (lowest_energy, index_k) )
		#lines = all_occurrences(self.outputfile, "Standard Nuclear Orientation")
		outf = open(self.outputfile, "r")
		atmsym = []
		coords = []
		for line in outf.readlines()[line1:line1+self.natoms]:
			columns = line.split()
			atmsym.append(columns[1])
			coords.append(float(columns[2]))
			coords.append(float(columns[3]))
			coords.append(float(columns[4]))
		if len(atmsym) == self.natoms and len(coords) == 3*self.natoms:
			self.lowest_energy_geometry = Molecule()
			self.lowest_energy_geometry.setup(atmsym, coords)
			self.has_last_geometry = 1


	def get_all_geometries(self, xyzfile='none'):
		lines = all_occurrences(self.outputfile, "Standard Nuclear Orientation")
		self.ngeom = len(lines)
		self.all_coords = []
		for k in range(0, self.ngeom):
			line1 = lines[k]
			outf = open(self.outputfile, "r")
			for line in outf.readlines()[line1+2:line1+self.natoms+2]:
				columns = line.split()
				self.all_coords.append(float(columns[2]))
				self.all_coords.append(float(columns[3]))
				self.all_coords.append(float(columns[4]))
			outf.close()
		self.all_energies = get_all_values(self.outputfile, "Total energy in the final basis set =", -1)
		if xyzfile != "none":
			xyzf = open(xyzfile, "w")
			for k in range(0, self.ngeom):
        			mol = Molecule()
        			mol.setup(qcout.input_geometry.atmnum, qcout.all_coords[k*3*natoms:(k+1)*3*natoms])
        			mol.mprint_with_energy(xyzfile, float(qcout.all_energies[k]))

	def get_scan_geometries(self, xyzfile='none'):
		all_energies = get_all_values(self.outputfile, "PES scan, value:", -1)
		ngeom = len(all_energies)
		all_coords = []
		outf = open(self.outputfile, "r")
		start = -1
		for line in outf.readlines():
			if "Coordinates (Angstroms)" in line:
				start = 1
				coords = []
			if start >= 3:
				columns = line.split()
				coords.append(float(columns[2]))
				coords.append(float(columns[3]))
				coords.append(float(columns[4]))
			if start == self.natoms + 2:
				start = -1
			if start >= 1:
				start += 1
			if "PES scan, value:" in line:
				for k in range(0, 3*self.natoms):
					all_coords.append(coords[k])
		outf.close()
		if xyzfile != "none":
			xyzf = open(xyzfile, "w")
			for k in range(0, ngeom):
        			mol = Molecule()
        			mol.setup(self.input_geometry.atmsym, all_coords[k*3*self.natoms:(k+1)*3*self.natoms])
        			mol.mprint_with_energy(xyzfile, float(all_energies[k]), 'append')

	def get_scf_energy(self):
		self.scf_energy = float(get_last_value(self.outputfile, "SCF   energy in the final basis set =", -1))



	def get_optimized_energy(self):
		self.optimized_energy = float(get_last_value(self.outputfile, "Final energy is", -1))

	def get_r12mr34_distance(self):
                self.r12mr34_distance = float(get_last_value(self.outputfile, "I:1 Type:2 Value", -4))



	def get_all_optimized_energies(self):
		energies = get_all_values(self.outputfile, "Final energy is", -1)
		self.all_optimized_energies = []
		for k in range(0, len(energies)):
			self.all_optimized_energies.append(float(energies[k]))


	def get_excited_state_energies(self):
		energies = get_all_values(self.outputfile, "Total energy for state   ", -1)
		self.excited_state_energies = []
		for k in range(0, len(energies)):
			self.excited_state_energies.append(float(energies[k]))


	def get_first_excited_state_energies(self):
		energies = get_all_values(self.outputfile, "Total energy for state   1", -1)
		self.first_excited_state_energies = []
		for k in range(0, len(energies)):
			self.first_excited_state_energies.append(float(energies[k]))


	def get_total_ground_state_energies(self):
		energies = get_all_values(self.outputfile, "Total energy in the final basis set", -1)
		self.total_ground_state_energies = []
		for k in range(0, len(energies)):
			self.total_ground_state_energies.append(float(energies[k]))


	def get_enthalpy_entropy(self):
		self.total_enthalpy = float(get_last_value(self.outputfile, "Total Enthalpy:", -2))
		self.total_entropy  = float(get_last_value(self.outputfile, "Total Entropy:", -2))

	def get_solvation_energy(self):
		self.solvation_energy  = float(get_last_value(self.outputfile, "(9) = (6) - (0)", -2))

	def get_grimme_d_dispersion_energy(self):
		self.grimme_d_dispersion_energy  = float(get_first_value(self.outputfile, "Empirical dispersion =", -2))

	def get_grimme_d3_dispersion_energy(self):
		self.grimme_d3_dispersion_energy  = float(get_first_value(self.outputfile, "D3 energy without 3body term =", -2))*AU2KCAL

	def get_grimme_d3bj_dispersion_energy(self):
		self.grimme_d3bj_dispersion_energy  = float(get_first_value(self.outputfile, "D3 energy without 3body term =", -2))*AU2KCAL

	def get_xdm_dispersion_energy(self):
		self.xdm_dispersion_energy  = float(get_last_value(self.outputfile, "Evdw(C6) =", -2))

	def get_mbd_dispersion_energy(self):
		self.mbd_dispersion_energy  = float(get_last_value(self.outputfile, "MBD Energy is", -1))*AU2KCAL

	def get_vv10_dispersion_energy(self):
		self.vv10_dispersion_energy  = float(get_last_value(self.outputfile, "Nonlocal correlation =", -1))*AU2KCAL

	def get_rimp2_dispersion_energy(self):
		energies = get_all_values(self.outputfile, "Total  RIMP2   correlation energy =", -2)
		print("energies=", energies)
		self.rimp2_dispersion_energy  = (float(energies[2])-float(energies[0])-float(energies[1]))*AU2KCAL

	def get_orbital_energies(self):
		line1 = last_occurrence(self.outputfile, "Orbital Energies (a.u.)")
		line2 = last_occurrence(self.outputfile, "Ground-State Mulliken")
		#print("line1: ", line1, "line2: ", line2)
		outf = open(self.outputfile, "r")
		self.alpha_orbital_energies = []
		for line in outf.readlines()[line1+4:line2-3]:
			#print("line=", line)
			cols = line.split()
			ncols = len(cols)
			if ncols > 0 and cols[0] != '--':
				for k in range(0, ncols):
					nch = len(cols[k])
					if nch > 8 or cols[k] == '*******':
						neles = int((nch+2)/9)
						for m in range(0, neles):
							self.alpha_orbital_energies.append(0.0)
					else:
						self.alpha_orbital_energies.append(float(cols[k]))

	def get_excitation_energies(self):
		energies = get_all_values(self.outputfile, "Excited state", -1)
		self.excitation_energies = []
		for k in range(0, len(energies)):
			self.excitation_energies.append(float(energies[k]))
		strengths = get_all_values(self.outputfile, "Strength   :", -1)
		self.oscillator_strength = []
		for k in range(0, len(strengths)):
			self.oscillator_strength.append(float(strengths[k]))


	def get_first_excitation_energies(self):
		energies = get_all_values(self.outputfile, "Excited state   1", -1)
		self.first_excitation_energies = []
		for k in range(0, len(energies)):
			self.first_excitation_energies.append(float(energies[k]))


	def get_ptss_excitation_energies(self):
		energies = get_all_values(self.outputfile, "Total  1st-order corrected excitation energy  =", -2)
		self.ptss_excitation_energies = []
		for k in range(0, len(energies)):
			self.ptss_excitation_energies.append(float(energies[k]))



	def get_excitation_energies_long(self):
		fields = get_all_values(self.outputfile, "Excited state", -5)
		fields2 = get_all_values(self.outputfile, "Excited state", -1)
		nstates = len(fields)
		self.excitation_energies = []
		self.oscillator_strength = []
		for k in range(0, nstates):
			self.excitation_energies.append(float(fields[k][:-2]))
			self.oscillator_strength.append(float(fields2[k][:-1]))

	def get_strengths(self):
		strengths = get_all_values(self.outputfile, "Strength",-1)
		self.strengths = []
		for k in range(0, len(strengths)):
			self.strengths.append(float(strengths[k]))


	def get_amplitudes(self):
		outf = open(self.outputfile, "r")
		line1 = last_occurrence(self.outputfile, "Excited state   1:")
		line2 = last_occurrence(self.outputfile, "Excited state   2:")
		self.nalpha = int(get_first_value_2(self.outputfile, "alpha", "beta", 2))
		print("nalpha=", int(get_first_value_2(self.outputfile, "alpha", "beta", 2)))
		self.amplitudes = []
		print("line1: ", line1, "line2: ", line2)
		for line in outf.readlines()[line1+5:line2]:
			columns = line.split()
			amplitudes = get_all_values(self.outputfile,"X: D( "+str(self.nalpha)+") --> V(  1)", -1)
			for k in range(0, len(amplitudes)):
				self.amplitudes.append(float(amplitudes[k]))

	def get_homo_lumo_gap(self):
		self.get_orbital_energies()
		return self.alpha_orbital_energies[self.nalpha]-self.alpha_orbital_energies[self.nalpha-1]

	def get_charges(self):
		line1 = last_occurrence(self.outputfile, "Ground-State Mulliken Net Atomic Charges")+3
		outf  = open(self.outputfile, "r")
		self.mulliken_charges = []
		for line in outf.readlines()[line1:line1+self.natoms]:
			#print("line=", line)
			cols = line.split()
			self.mulliken_charges.append(float(cols[2]))
		#print("mullken.charges = ", self.mulliken_charges)
		outf.close()
		lines = []
		if last_occurrence(self.outputfile, "TDA Excited State") > 0:
			lines = all_occurrences(self.outputfile, "TDA Excited State")
		elif last_occurrence(self.outputfile, "RPA Excited State") > 0:
			lines = all_occurrences(self.outputfile, "RPA Excited State")
		#lines = all_occurrences(self.outputfile, "TDA Excited State")
		#print("lines: ", lines)
		for k in range(0, len(lines)):
			line1 = int(lines[k])+3
			#print("k: ", k, "line1:", line1)
			outf  = open(self.outputfile, "r")
			for line in outf.readlines()[line1:line1+self.natoms]:
				#print("line:", line)
				cols = line.split()
				self.mulliken_charges.append(float(cols[2]))
			outf.close()
		#matrix_print_1d(self.mulliken_charges, self.natoms, 6, 6, "mulliken charges")
		line1 = last_occurrence(self.outputfile, "Merz-Kollman ESP Net Atomic Charges")+3
		if line1 > 30:
			outf = open(self.outputfile, "r")
			#print("line1: ", line1)
			self.esp_charges = []
			for line in outf.readlines()[line1:line1+self.natoms]:
				#print("line=", line)
				cols = line.split()
				self.esp_charges.append(float(cols[2]))
			outf.close()
		line2 = last_occurrence(self.outputfile, "ESP charges for excited states")+1
		if line2 > 30:
			outf = open(self.outputfile, "r")
			nstates = int(len(self.mulliken_charges)/self.natoms) - 1
			if nstates == 0:
				nstates = int(get_last_value(self.outputfile, "NStates=", -1)) - 1
			#print("nstates=", nstates)
			charges = self.get_real_block(self.natoms, nstates, "ESP charges for excited states")
			for j in range(0, nstates):
				for i in range(0, self.natoms):
					#self.esp_charges.append(charges[j+i*nstates])
					self.esp_charges.append(charges[j,i])
			#matrix_print_1d(self.esp_charges, self.natoms, nstates+1, 6, "esp charges")

	# dipole moment
	def get_dipole(self):
		values = np.zeros(4)
		line1 = last_occurrence(self.outputfile, "Dipole Moment")
		outf = open(self.outputfile, "r")
		for line in outf.readlines()[line1:line1+2]:
			#print('line=',line)
			cols = line.split()
			if cols[0] == 'X':
				values[0] = float(cols[1])
				values[1] = float(cols[3])
				values[2] = float(cols[5])
			elif cols[0] == 'Tot':
				values[3] = float(cols[1])
		outf.close()
		return values

	# m is the leading dimension of the array
	# should remove this later, and use the generic one below
	def get_real_square_block(self, array, M, Words):
		words = Words.split()
		len_words = len(words)

		outf = open(self.outputfile, "r")
		start = 0
		column_indices =[0, 0, 0, 0, 0, 0]
		for line in outf.readlines():
			columns= line.split()
			ncols = len(columns)
			if start == 1 and ncols >= 1:
				#print('line=', line, is_integer(columns[0]), is_integer(columns[1]))
				if is_integer(columns[0]) == 0:    #end of block
					start = 2
				elif ncols > 1 and is_integer(columns[1]) == 0: # real data
					row_index = int(columns[0]) - 1
					#print('row_index=', row_index)
					for k in range(1, ncols):
						array[row_index + column_indices[k-1] * M] = float(columns[k])
				else:
					for k in range(0, ncols):
						column_indices[k] = int(columns[k])-1
					#print('column_indices=', column_indices)
			if start == 0 and ncols == len_words:
				is_block = 1
				for k in range(0, ncols):
					if columns[k] != words[k]: is_block = 0
				if is_block == 1:
					start = 1

	# read a block of data printed using MatPrint
	def get_real_block(self, M, N, Words):
		words = Words.split()
		len_words = len(words)

		#print("N=", N, "M=", M)
		array = np.zeros((N, M))

		outf = open(self.outputfile, "r")
		start = 0
		column_indices =[0, 0, 0, 0, 0, 0]
		for line in outf.readlines():
			columns= line.split()
			ncols = len(columns)
			if start == 1 and ncols >= 1:
				#print('line=', line)
				if is_integer(columns[0]) == 0:    #end of block
					start = 2
				elif ncols > 1 and is_integer(columns[1]) == 0: # real data
					row_index = int(columns[0]) - 1
					#print('row_index=', row_index)
					for k in range(1, ncols):
						array[column_indices[k-1], row_index] = float(columns[k])
				else:
					for k in range(0, ncols):
						column_indices[k] = int(columns[k])-1
					#print('column_indices=', column_indices)
			if start == 0 and ncols == len_words:
				is_block = 1
				for k in range(0, ncols):
					if columns[k] != words[k]: is_block = 0
				if is_block == 1:
					start = 1
		return array

	def get_z_matrix(self):
		line1 = last_occurrence(self.outputfile, "Z-matrix Print:")
		outf = open(self.outputfile, "r")
		self.z_elements = []
		self.z_values = np.zeros(6*self.natoms)
		icount = 0
		for line in outf.readlines()[line1+2:line1+2+self.natoms]:
			cols = line.split()
			ncols = len(cols)
			self.z_elements.append(cols[0])
			for k in range(1, ncols-1, 2):
				self.z_values[6*icount+k-1] = float(cols[k])+0.001
				self.z_values[6*icount+k]   = float(cols[k+1])
			icount += 1

		#print(self.z_elements)
		#matrix_print_1d(self.z_values, 6, self.natoms, 6, "z values")

	def find_harmonic_frequencies(self):
		need_words = "Frequency:"
		loop_index = 0
		self.harmonic_frequencies = []
		with open(self.outputfile, 'r', encoding="ISO-8859-1") as infile:
			for line in infile:
				if line.find(need_words) >= 0:
					catch_words = line.split()
					for k in range(0, len(catch_words)-1):
						if catch_words[k+1] == "********": catch_words[k+1] = 0.0
						self.harmonic_frequencies.append(float(catch_words[k+1]))
						loop_index += 1
		print("frequencies:", self.harmonic_frequencies)

	def find_intensities(self, ir_or_raman):
		if ir_or_raman == 'IR':
			need_words = 'IR Intens:'
			self.ir_intensities = []
		else:
			need_words = 'Raman Intens:'
			self.raman_intensities = []
		loop_index = 0
		with open(self.outputfile, 'r', encoding="ISO-8859-1") as infile:
			for line in infile:
				if line.find(need_words) >= 0:
					catch_words = line.split()
					for k in range(0, len(catch_words)-2):
						if ir_or_raman == 'IR':
							self.ir_intensities.append(float(catch_words[k+2]))
						else:
							self.raman_intensities.append(float(catch_words[k+2]))
					loop_index += 1
		if ir_or_raman == 'IR':
			print("intensities:", self.ir_intensities)
		else:
			print("intensities:", self.raman_intensities)


# -*- coding: utf-8 -*-
#
# Copyright (c), the QMAnalysis development team
#
#
#

""" Common functions for making plot """

import numpy, sys, os
from .plot_util import *
from .util import *
from .qchem_out_parser import *
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

def make_jmol_isosurface_plot(img_file, cube_file, cutoff=0.02):

	spt = open("cube.spt", "w")
	spt.write("background = white\n")
	spt.write("load "+cube_file+"\n")
	#spt.write("rotate x 45\n")
	#spt.write("rotate x -30\n")
	#spt.write("rotate z 60\n")
	#spt.write("rotate y 100\n")
	if os.path.isfile("cube0.spt"):
		tfile = open("cube0.spt", "r")
		for iline in tfile.readlines():
			spt.write(iline)
		tfile.close()
	spt.write("isosurface cutoff "+str(cutoff)+" sign color translucent "+cube_file+"\n")
	spt.write("write image "+img_file+"\n")
	spt.close()
	os.system("./jmol -n cube.spt")

def maximum_degeneracy(energies):
	norb = len(energies)
	maxd = 1
	for k in range(0, norb):
		ndegeneracy = 1
		m = k+1
		while(m < norb and energies[m] - energies[k] < 0.0001):
			ndegeneracy += 1
			m += 1
		if ndegeneracy > maxd: maxd = ndegeneracy
	return maxd

def make_orbital_energy_plot(plotfile, energies, occupancies):

	norb   = len(energies)
	emin = np.amin(energies)
	emax = np.amax(energies)
	maxdegeneracy = maximum_degeneracy(energies)
	#print("maxdegeneracy:", maxdegeneracy)

	width  = 30
	total_width = width * (maxdegeneracy + 0.2)
	height = 300
	fig = plt.figure(figsize=(2*(maxdegeneracy+0.2),6), dpi=600)

	text_skip = []
	for k in range(0, norb):
		text_skip.append(0)
	xtext = ((maxdegeneracy - 1) * 0.5 + 0.35) * width + 2

	k = 0
	while k <= norb-1:

		#print("k:", k)
		ndegeneracy = 1
		m = k+1
		while(m < norb and energies[m] - energies[k] < 0.0001):
			ndegeneracy += 1
			m += 1
		for m in range(0, ndegeneracy):
			xmin  = (m - (ndegeneracy-1) * 0.5 - 0.22) * width
			xmax  = (m - (ndegeneracy-1) * 0.5 + 0.22) * width
			x1    = (m - (ndegeneracy-1) * 0.5 - 0.05) * width  #arrow 1
			x2    = (m - (ndegeneracy-1) * 0.5 + 0.05) * width  #arrow 2
			y    = (energies[k] - emin) / (emax - emin) * height
			plt.plot([xmin, xmax], [y, y], color='black')
			if occupancies[k+m] == 1:
				add_arrow(fig, x1, y-3, 0,  6, 'red', 'none')
			elif occupancies[k+m] == 2:
				add_arrow(fig, x1, y-3, 0,  6, 'blue', 'none')
				add_arrow(fig, x2, y+3, 0, -6, 'blue', 'none')
			#print("nd:", ndegeneracy, "m:", m, "xmin, xmax=", xmin, xmax)
		if text_skip[k] == 0:
			#print("energy: ", energies[k], float_str(energies[k], 3))
			if energies[k] < 0:
				if (energies[k+ndegeneracy] - energies[k]) / (emax-emin) > 0.015:
					plt.text(xtext, y-2, float_str(energies[k]), fontsize=8, color='blue')
				else:
					plt.text(xtext, y-2, float_str(energies[k])+', '+float_str(energies[k+ndegeneracy]), fontsize=8, color='blue')
					text_skip[k+ndegeneracy] = 1
			else:
				if k+ndegeneracy == norb or (energies[k+ndegeneracy] - energies[k]) / (emax-emin) > 0.015:
					plt.text(xtext, y-2, float_str(energies[k]), fontsize=8, color='red')
				else:
					plt.text(xtext, y-2, float_str(energies[k])+', '+float_str(energies[k+ndegeneracy]), fontsize=8, color='red')
					text_skip[k+ndegeneracy] = 1
		k += ndegeneracy

	plt.xlim([-width*maxdegeneracy/2, width*(maxdegeneracy/2+0.6)])
	plt.axis('off')
	plt.tight_layout()
	plt.savefig(plotfile)

def fit_val(positions, heighs, broaden):
	num_vib = len(positions)
	x_min = positions[0] - 50
	if x_min < 0:
		x_min = 0
	x_max = positions[num_vib-1] + 50
	if x_max > 4000:
		x_max = 4000
	ix = np.linspace(x_min, x_max, int(x_max-x_min))

	iy = 0.0
	for peak in range(0, num_vib):
		iy += 2.51225 * broaden * heighs[peak] \
			 * mlab.normpdf(ix,positions[peak],broaden)
	return (ix, iy)

def plot_ir_or_raman_spectra(plotfile, qchemoutfile, ir_or_raman, broaden):
	qcout = QChemOut(qchemoutfile)
	qcout.find_harmonic_frequencies()
	qcout.find_intensities('Raman')

	mode = qcout.natoms * 3 - 6
	peak_centers = qcout.harmonic_frequencies
	peak_intens = qcout.raman_intensities
	#find_intens(peak_intens, "IR Intens:", out_file)
	#find_intens(peak_intens, "Raman Intens:", out_file)
	ix, iy = fit_val(peak_centers, peak_intens, broaden)
	plt.plot(ix, iy)
	plt.xlim(800,3600)
	plt.ylim(0, 1200)
	plt.xticks([1000,1500,2000,2500,3000,3500], size=14)
	plt.yticks([0,400,800,1200], size=14)
	plt.xlabel("Frequency (cm$^{-1}$)",fontsize=16)
	plt.ylabel("Intensity",fontsize=16)
	plt.tight_layout()
	plt.savefig(plotfile)

#This adds a subplot
def plot_ir_or_raman_spectra2(ax, qchemoutfile, ir_or_raman, broaden):
	qcout = QChemOut(qchemoutfile)
	qcout.find_harmonic_frequencies()
	qcout.find_intensities('Raman')

	mode = qcout.natoms * 3 - 6
	peak_centers = qcout.harmonic_frequencies
	peak_intens = qcout.raman_intensities
	#find_intens(peak_intens, "IR Intens:", out_file)
	#find_intens(peak_intens, "Raman Intens:", out_file)
	ix, iy = fit_val(peak_centers, peak_intens, broaden)
	ax.plot(ix, iy, color='red')
	ax.set_xlim(3700, 800)
	ax.set_ylim(0, 1200)
	ax.set_xticks([1000,1500,2000,2500,3000,3500])
	ax.set_yticks([0,400,800,1200])
	ax.tick_params(axis="x", labelsize=7)
	ax.tick_params(axis="y", labelsize=7)
	#ax.set_xlabel("Frequency (cm$^{-1}$)",fontsize=16)
	#ax.set_ylabel("Intensity",fontsize=16)

from .molecule import Molecule
from .external_charges import *
from .util import *
import sys, os, string

class QChemRem():

	""" set up QChem rem variables """
	def __init__(self, dict=None):
		self.dict = {'jobtype' : 'sp', 
			     'method' : 'b3lyp', 
			     'thresh' : '14',
			     'scf_convergence' : '8',
			     'mem_static' : '500',
			     'xc_grid' : '1',
			     'basis' : '6-31G*'}
		for key in dict:
			self.dict[key] = dict[key]
		if 'cis_n_roots' in dict.keys():
			self.dict['set_iter'] = '100'


def qchem_input_generator(inpfile, molecule, rem, append='none'):

	if append == 'none':
		inpf = open(inpfile, "w")
	elif append == 'append' or append == 'read-molecule':
		inpf = open(inpfile, "a")
		inpf.write("@@@ \n\n")
	inpf.write("$molecule\n")
	if append != 'read-molecule':
		inpf.write("%d %d\n" % (molecule.charge, molecule.multiplicity))
		for k in range(0, molecule.natoms):	
			inpf.write(" %3s" %(molecule.atmsym[k]))
			for m in range(0, 3):
				inpf.write(" %15.10f" % (molecule.coords[3*k+m]*BOHR))
			inpf.write("\n")
	else: 
		inpf.write("read\n")
	inpf.write("$end\n\n")

	inpf.write("$rem\n")
	for key in sorted(rem.dict):
		inpf.write("%s %s\n" % (key, rem.dict[key]))
	inpf.write("$end\n\n")
	inpf.close()

def qchem_input_generator_with_external_charges(inpfile, molecule, rem, extCharges, potential='none'):
	qchem_input_generator(inpfile, molecule, rem)
	if extCharges.nchgs > 0:
		append_external_charges(inpfile, extCharges)
	if potential != 'none':
		append_atom_site_potential(inpfile, potential)

def qchem_input_generator_selected_atoms(inpfile, molecule, rem, list, append='none'):

	if append == 'none':
		inpf = open(inpfile, "w")
	elif append == 'append' or append == 'read-molecule':
		inpf = open(inpfile, "a")
		inpf.write("@@@ \n\n")
	inpf.write("$molecule\n")
	inpf.write("%d %d\n" % (molecule.charge, molecule.multiplicity))
	for k in range(0, molecule.natoms):	
		if list[k] == 1:
			inpf.write(" %3s" %(molecule.atmsym[k]))
			for m in range(0, 3):
				inpf.write(" %15.10f" % (molecule.coords[3*k+m]*BOHR))
			inpf.write("\n")
	inpf.write("$end\n\n")

	inpf.write("$rem\n")
	for key in sorted(rem.dict):
		inpf.write("%s %s\n" % (key, rem.dict[key]))
	inpf.write("$end\n\n")
	inpf.close()

#two jobs in the same input
def qchem_2_input_generator(inpfile, molecule, rem0, rem):

	inpf = open(inpfile, "w")
	inpf.write("$molecule\n")
	inpf.write("%d %d\n" % (molecule.charge, molecule.multiplicity))
	for k in range(0, molecule.natoms):	
		inpf.write(" %3s" %(molecule.atmsym[k]))
		for m in range(0, 3):
			inpf.write(" %15.10f" % (molecule.coords[3*k+m]*BOHR))
		inpf.write("\n")
	inpf.write("$end\n\n")

	inpf.write("$rem\n")
	for key in sorted(rem0.dict):
		inpf.write("%s %s\n" % (key, rem0.dict[key]))
	inpf.write("$end\n\n")

	inpf.write("@@@\n\n")
	inpf.write("$molecule\n")
	inpf.write("READ\n")
	inpf.write("$end\n\n")

	inpf.write("$rem\n")
	for key in sorted(rem.dict):
		inpf.write("%s %s\n" % (key, rem.dict[key]))
	inpf.write("$end\n\n")

	inpf.close()

#two jobs in the same input with external charges
def qchem_2_input_generator_with_external_charges(inpfile, molecule, rem0, rem, extCharges0, extCharges, potential):

	inpf = open(inpfile, "w")
	inpf.write("$molecule\n")
	inpf.write("%d %d\n" % (molecule.charge, molecule.multiplicity))
	for k in range(0, molecule.natoms):	
		inpf.write(" %3s" %(molecule.atmsym[k]))
		for m in range(0, 3):
			inpf.write(" %15.10f" % (molecule.coords[3*k+m]*BOHR))
		inpf.write("\n")
	inpf.write("$end\n\n")

	inpf.write("$rem\n")
	for key in sorted(rem0.dict):
		inpf.write("%s %s\n" % (key, rem0.dict[key]))
	inpf.write("$end\n\n")

	inpf.write("$external_charges\n")
	for k in range(0, extCharges0.nchgs):	
		for m in range(0, 3):
			inpf.write(" %15.10f" % (extCharges0.coords[3*k+m]*BOHR))
		inpf.write(" %15.10f\n" % (extCharges0.cutoff_charges[k]))
	inpf.write("$end\n\n")

	if len(potential) > 0: 
		inpf.write("$atom_site_potential\n")
		for k in range(0, int(len(potential)/4)):
			inpf.write("%3d %12.7f " % (k+1, potential[4*k+0]))
			for m in range(0, 3):
				inpf.write(" %12.7f " % (potential[4*k+m+1]))
			inpf.write("\n")
		inpf.write("$end\n\n")

	inpf.write("@@@\n\n")
	inpf.write("$molecule\n")
	inpf.write("READ\n")
	inpf.write("$end\n\n")

	inpf.write("$rem\n")
	for key in sorted(rem.dict):
		inpf.write("%s %s\n" % (key, rem.dict[key]))
	inpf.write("$end\n\n")

	inpf.write("$external_charges\n")
	for k in range(0, extCharges.nchgs):	
		for m in range(0, 3):
			inpf.write(" %15.10f" % (extCharges.coords[3*k+m]*BOHR))
		inpf.write(" %15.10f\n" % (extCharges.cutoff_charges[k]))
	inpf.write("$end\n\n")

	inpf.close()

def qchem_input_generator_with_remfile(inpfile, molecule, remfile):

	inpf = open(inpfile, "w")
	inpf.write("$molecule\n")
	inpf.write("%d %d\n" % (molecule.charge, molecule.multiplicity))
	for k in range(0, molecule.natoms):	
		inpf.write(" %3s" %(molecule.atmsym[k]))
		for m in range(0, 3):
			inpf.write(" %15.10f" % (molecule.coords[3*k+m]*BOHR))
		inpf.write("\n")
	inpf.write("$end\n\n")

	remf = open(remfile, "r")
	for line in remf.readlines():
		inpf.write("%s" % (line))
	remf.close()

	inpf.close()

#  molecular information from z matrix 
def qchem_input_generator_with_zmatrix(inpfile, molecule, rem, append='none'):

	if append == 'none':
		inpf = open(inpfile, "w")
	elif append == 'append' or append == 'read-molecule':
		inpf = open(inpfile, "a")
		inpf.write("@@@ \n\n")
	inpf.write("$molecule\n")
	if append != 'read-molecule':
		inpf.write("%d %d\n" % (molecule.charge, molecule.multiplicity))
		for k in range(0, molecule.natoms):	
			inpf.write(" %3s" %(molecule.z_elements[k]))
			for m in range(0, 3):
				if molecule.z_values[6*k+2*m] > 0: 
					inpf.write(" %3d %15.10f" % (molecule.z_values[6*k+2*m], molecule.z_values[6*k+2*m+1]))
			inpf.write("\n")
	else:
		inpf.write("read\n")
	inpf.write("$end\n\n")

	inpf.write("$rem\n")
	for key in sorted(rem.dict):
		inpf.write("%s %s\n" % (key, rem.dict[key]))
	inpf.write("$end\n\n")

	inpf.close()

#bsse style
def qchem_frag_input_generator(inpfile, molecule, rem, append=False):

	if append:
		inpf = open(inpfile, "a")
		inpf.write("\n@@@\n\n")
	else:
		inpf = open(inpfile, "w")
	inpf.write("$molecule\n")
	inpf.write("%d %d\n" % (molecule.charge, molecule.multiplicity))
	for ifrag in range(0, molecule.nfrag):
		inpf.write("--\n")
		if not hasattr(molecule, 'frag_charges'):
			inpf.write("0 1\n")
		else:
			inpf.write("%d 1\n" % molecule.frag_charges[ifrag])
		for k in range(0, molecule.natoms):	
			if molecule.frag_index[k] == ifrag:
				inpf.write(" %3s" %(molecule.atmsym[k]))
				for m in range(0, 3):
					inpf.write(" %15.10f" % (molecule.coords[3*k+m]*BOHR))
				inpf.write("\n")
	inpf.write("$end\n\n")

	inpf.write("$rem\n")
	for key in sorted(rem.dict):
		inpf.write("%s %s\n" % (key, rem.dict[key]))
	inpf.write("$end\n\n")
	inpf.close()


#frag1, frag2, ..., total complex
def qchem_frag_seq_input_generator(inpfile, molecule, rem):

	inpf = open(inpfile, "w")
	for ifrag in range(0, molecule.nfrag+1):
		inpf.write("$molecule\n")
		inpf.write("0 1\n")
		for k in range(0, molecule.natoms):	
			if molecule.frag_index[k] == ifrag or ifrag == molecule.nfrag:
				inpf.write(" %3s" %(molecule.atmsym[k]))
				for m in range(0, 3):
					inpf.write(" %15.10f" % (molecule.coords[3*k+m]*BOHR))
				inpf.write("\n")
		inpf.write("$end\n\n")

		inpf.write("$rem\n")
		for key in sorted(rem.dict):
			inpf.write("%s %s\n" % (key, rem.dict[key]))
		inpf.write("$end\n\n")
		if ifrag != molecule.nfrag: 
			inpf.write("@@@\n\n")
	inpf.close()

def append_classical_water(inpfile, molecule, iatom_first_water, model):
	inpf = open(inpfile, "a")
	if model == 'tip3p' or model == 'tip5p':
		inpf.write("$external_charges\n")
	#efp2 uses HF/6-31G* parameters for TIP3P-geometry water molecules
	#efp3 uses HF/6-31+G* parameters for TIP3P-geometry water molecules
	elif model == 'efp' or model == 'efp2' or model == 'efp3' or model == 'amoeba':
		inpf.write("$efp_fragments\n")
	n_water = int((molecule.natoms - iatom_first_water+1)/3)
	for m in range(0, n_water):
		if model == 'efp':
			inpf.write("WATER_L\n")
		elif model == 'efp2':
			inpf.write("TIP3P_L\n")
		elif model == 'efp3':
			inpf.write("TIP3PD_L\n")
		elif model == 'amoeba':
			inpf.write("AMOEBAWATER_L\n")
		if model != 'tip5p':
			for i in range(0, 3):  # three atoms
				iatom = iatom_first_water-1+3*m+i
				#sanity check
				if i == 0 and molecule.atmnum[iatom] != 8:
					inpf.write("not an oxygen")
				elif i > 0 and molecule.atmnum[iatom] != 1:
					inpf.write("not an hydrogen")
				if model == 'efp' or model == 'efp2' or model == 'efp3' or model == 'amoeba':
					if i == 0: inpf.write("O ")
					else: inpf.write("H ")
				for j in range(0, 3): 
					inpf.write("%15.8f " % (molecule.coords[3*iatom+j]*BOHR))
				if model == 'tip3p': 
					if i == 0: inpf.write("    -0.834\n")
					else: inpf.write("     0.417\n")
				else: 
					inpf.write("\n")
		else:
			#tip5p
			atom1 = iatom_first_water+3*m
			charges, xyz = tip5p_data(molecule.coords, atom1)
			for i in range(0, 4): 
				for j in range(0, 3): 
					inpf.write("%15.8f " % (xyz[j+i*3]*BOHR))
				inpf.write("%15.3f \n" % (charges[i]))
	inpf.write("$end\n\n")

#print charges within a cutoff (in Angstrom)
def append_external_charges(inpfile, extCharges, cutoff=10000.0):
	inpf = open(inpfile, "a")
	inpf.write("$external_charges\n")
	for k in range(0, extCharges.nchgs):	
		if extCharges.distances[k] < cutoff:
			for m in range(0, 3):
				inpf.write(" %15.10f" % (extCharges.coords[3*k+m]*BOHR))
			inpf.write(" %15.10f\n" % (extCharges.cutoff_charges[k]))
	inpf.write("$end\n\n")

def append_atom_site_potential(inpfile, potential):
	inpf = open(inpfile, "a")
	inpf.write("$atom_site_potential\n")
	for k in range(0, int(len(potential)/4)):
		inpf.write("%3d %12.7f " % (k+1, potential[4*k+0]))
		for m in range(0, 3):
			inpf.write(" %12.7f " % (potential[4*k+m+1]))
		inpf.write("\n")
	inpf.write("$end\n\n")

def append_geometry_constraints(inpfile, constraints):
	inpf = open(inpfile, "a")
	inpf.write("$opt\n")
	inpf.write("CONSTRAINT\n")
	print("constraints=", constraints)
	for k in range(0, len(constraints)):
		if constraints[k][0] == 'stre':
			inpf.write("%4s %3d %3d %12.7f\n" % (constraints[k][0], constraints[k][1],
					 constraints[k][2], constraints[k][3]))
		elif constraints[k][0] == 'bend':
			inpf.write("%4s %3d %3d %3d %12.7f\n" % (constraints[k][0], constraints[k][1],
					 constraints[k][2], constraints[k][3], constraints[k][4]))
		elif constraints[k][0] == 'tors':
			inpf.write("%4s %3d %3d %3d %3d %12.7f\n" % (constraints[k][0], constraints[k][1],
					 constraints[k][2], constraints[k][3], constraints[k][4], constraints[k][5]))
	inpf.write("ENDCONSTRAINT\n")
	inpf.write("$end\n\n")

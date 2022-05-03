from .molecule import Molecule
from .external_charges import *
from .util import *
import sys
write = sys.stdout.write

def input_section(filename, section_name):
	first_line = first_occurrence(filename, section_name)
	#print("first_line=", first_line)
	all_lines = all_occurrences(filename, "\$end")
	#print("all_lines: ", all_lines)
	for k in range(0, len(all_lines)):
		if all_lines[k] > first_line:
			last_line = all_lines[k]
			break
	#print("first_line: ", first_line, "last_line", last_line)
	return first_line, last_line

def read_molecule(filename):

	first_line, last_line = input_section(filename, "\$molecule")
	inpf = open(filename, "r")
	atomsym = []
	xyz     = []
	for line in inpf.readlines()[first_line-1:last_line]:
		columns = line.split()
		ncols = len(columns)
		if ncols == 4:
			atomsym.append(columns[0])
			xyz.append(float(columns[1]))
			xyz.append(float(columns[2]))
			xyz.append(float(columns[3]))

	molecule = Molecule()
	molecule.setup(atomsym, xyz)

	return molecule

def read_rem(filename):
	first_line, last_line = input_section(filename, "\$rem")
	inpf = open(filename, "r")
	options = {}
	for line in inpf.readlines()[first_line:last_line-1]:
		#print('line=', line)	
		columns = line.split()
		ncols = len(columns)
		if ncols == 2:
			options[columns[0].lower()] = columns[1].lower()

	return options

def read_charges(filename):
	first_line, last_line = input_section(filename, "\$external_charges")
	inpf = open(filename, "r")
	charges = []
	xyz     = []
	icount = 0
	for line in inpf.readlines()[first_line-1:last_line]:
		columns = line.split()
		if int(icount/1000000) * 1000000 == icount: print("icount=", icount)
		ncols = len(columns)
		if ncols == 4:
			xyz.append(float(columns[0]))
			xyz.append(float(columns[1]))
			xyz.append(float(columns[2]))
			charges.append(float(columns[3]))
		icount += 1
	print("number of external charges: ", len(charges))

	extCharges = ExternalCharges()
	extCharges.setup(xyz, charges)
	return extCharges

def read_qchem_input(filename):

	molecule = read_molecule(filename)
	options = read_rem(filename)
	extCharges = read_charges(filename)

	return molecule, extCharges, options

	

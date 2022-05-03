from .molecule import Molecule
from .util import *
import sys, os, string

class GaussianRoute():

	""" set up Gaussian rem variables """
	def __init__(self, dict=None):
		#default settings
		self.dict = {'functional' : 'b3lyp', 
			     'basis' : '6-31+G**'}
		for key in dict:
			self.dict[key] = dict[key]
			if key == 'functional' and dict[key] == 'b3lyp' and ('empiricaldispersion' not in key):
				self.dict['empiricaldispersion']='GD3BJ'

def gaussian_input_generator(inpfile, molecule, route):


#%nprocshared=12
#%mem=12GB
#%chk=oxyL_SNSN_keto-O_exp_opt.chk
##p mPW1PW91 6-31++G(d,p)
#Int=UltraFineGrid
#SCF=VeryTight
#Opt=(VTight,MaxCycle=1024,MaxStep=5)
#SCRF=(SMD,Solvent=Water,Read)
#NoSymm
#Density=Current
#Pop=CHelpG
#

	inpf = open(inpfile, "w")

	if "nproc" in route.dict:
		inpf.write("%%nprocshared=%s\n" % route.dict["nproc"])
	else:
		inpf.write("%%nprocshared=1\n" % route.dict["nproc"])
	if "mem" in route.dict:
		inpf.write("%%mem=%s\n" % route.dict["mem"])
	else:
		inpf.write("%%mem=10MB\n" % route.dict["mem"])
	inpf.write("%%chk=%s\n" % (inpfile[:-3]+'chk'))
	inpf.write("#p %s %s\n" % (route.dict['functional'], route.dict['basis']))
	for key in sorted(route.dict):
		if key == 'sym':
			inpf.write("%s\n" % (route.dict[key]))
		elif key != 'nproc' and key != 'mem' and key != 'functional' and key != 'basis':
			inpf.write("%s=%s\n" % (key, route.dict[key]))

	inpf.write("\ntitle\n\n")

	inpf.write("%d %d \n" % (molecule.charge, molecule.multiplicity))
	for k in range(0, molecule.natoms):	
		inpf.write(" %3s" %(molecule.atmsym[k]))
		for m in range(0, 3):
			inpf.write(" %15.10f" % (molecule.coords[3*k+m]*BOHR))
		inpf.write("\n")
	inpf.write("\n")



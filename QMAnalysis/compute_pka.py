import sys, string
from qchem_out_parser import *
from util import *

def compute_pKa(outputfile1, outputfile2):

	""" outputfile1 must be the acid (HX)
	    outputfiel2 must be the base (X) """

	out1 = QChemOut()
	out1.setup(outputfile1)
	out1.get_solvation_energy(outputfile1)
	out1.get_enthalpy_entropy(outputfile1)

	out2 = QChemOut()
	out2.setup(outputfile2)
	out2.get_solvation_energy(outputfile2)
	out2.get_enthalpy_entropy(outputfile2)

	print 'acid data from', outputfile1
	print 'base data from', outputfile2
	print ' '

        G_proton = -4.39 
	dG_gas = (out2.optimized_energy - out1.optimized_energy) * AU2KCAL + G_proton
        dG_solv1 = out1.solvation_energy
        dG_solv2 = out2.solvation_energy
	dG_solvp = -265.9 

	dG = dG_gas + dG_solv2 + dG_solvp - dG_solv1
	pKa = dG * CAL2J * 1000 / ( 8.314 * 298.15) / 2.303

	print("acid,   energy: %12.7f au solvation: %9.4f kcal/mol" % (out1.optimized_energy, out1.solvation_energy))
	print("base,   energy: %12.7f au solvation: %9.4f kcal/mol" % (out2.optimized_energy, out2.solvation_energy))
	print("proton, energy: %6.3f kcal/mol solvation: %9.4f kcal/mol (reference values)" % (G_proton, dG_solvp))
	print("deltaG:         %6.3f kcal/mol (without enthalpy/entropy corrections) " % (dG))
	print("pKa:            %6.3f" % (pKa))
	print(" ")

	dG = dG + out2.total_enthalpy - out1.total_enthalpy - 298.15 * 0.001 * (out2.total_entropy - out1.total_entropy)
	pKa = dG * CAL2J * 1000 / ( 8.314 * 298.15) / 2.303

	print("acid, enthalpy: %6.3f kcal/mol   entropy: %9.3f kcal/mol/K " % (out1.total_enthalpy, out1.total_entropy))
	print("base, enthalpy: %6.3f kcal/mol   entropy: %9.3f kcal/mol/K " % (out2.total_enthalpy, out2.total_entropy))
	print("deltaG:         %6.3f kcal/mol (with enthalpy/entropy corrections)" % (dG))
	print("pKa w/ corr:    %6.3f" % (pKa))

	return pKa

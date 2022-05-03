import numpy as np
from .math_util import *

def get_dm_from_occMO(Coa, Cob):

        DMa = np.dot(Coa.transpose(), Coa)
        DMb = np.dot(Cob.transpose(), Cob)
        DM = np.add(DMa, DMb)

        #matrix_print_2d(DM, 6, "DM")
        return DM

def get_dm_from_nonorthogonal_occMO(Coa, Cob, overlap_matrix):

	Cmat = Coa
	nbas = Coa.shape[1]
	print("nbas: ", nbas)
	Smat = np.reshape(overlap_matrix, (-1, nbas))
	
	SC   = np.dot(Cmat, Smat)
	CtSC = np.dot(SC, Cmat.transpose())

	CtSC_invsqrt = symmetric_matrix_power(CtSC, -0.5)
	Coa_new = np.dot(CtSC_invsqrt, Cmat)
	DMa = np.dot(Coa_new.transpose(), Coa_new)
	Dmb = DMa
	DM = np.add(DMa, Dmb)
	return DM
	
def check_orbital_orthogonality(C, S, nbas, norb):

	Cmat = np.reshape(C, (-1, nbas))
	Smat = np.reshape(S, (-1, nbas))
	#matrix_print_2d(Cmat, 6, "Cmat")
	#matrix_print_2d(Smat, 6, "Smat")

	SC   = np.dot(Cmat, Smat)
	CtSC = np.dot(SC, Cmat.transpose())
	matrix_print_2d(CtSC, 10, "CtSC")

	#print("Cmat: shape:", Cmat.shape)
	#print("Smat: shape:", Smat.shape)
	#print("SC: shape:", SC.shape)
	#print("CtSC: shape:", CtSC.shape)

def check_density_matrix_idempotency(P, S, nbas):
	Pmat = np.reshape(P, (-1, nbas))
	Smat = np.reshape(S, (-1, nbas))
	
	PS  = np.dot(Smat, Pmat)
	PSP = np.dot(Pmat, PS)
	I   = PSP - Pmat
	#matrix_print_2d(I, 10, "PSP - P")
	vmax = 0.0
	for j in range(0, nbas):
		for i in range(0, nbas):
			if abs(I[i,j]) > vmax: 
				vmax = abs(I[i,j])
	print("Largest Element in PSP-P:", vmax)
			

def diagonalize_density_matrix(nbas, P, overlap_matrix):

	#matrix_print_2d(P, 6, "P")
	S = np.reshape(overlap_matrix, (-1, nbas))	
	Shalf  = linalg.sqrtm(S)
	SInvHalf = symmetric_matrix_power(S, -0.5)
	#Shalf = symmetric_matrix_sqrt(S)
	#matrix_print_2d(Shalf, 6, "Shalf")

	ShalfP = np.dot(P, Shalf)
	ShalfPShalf = np.dot(Shalf, ShalfP)
	#matrix_print_2d(ShalfPShalf, 6, "ShalfPShalf")

	eval, evec = linalg.eigh(ShalfPShalf)
	matrix_print_1d(eval, 1, nbas, 6, "eval")
	matrix_print_2d(evec, 6, "evec")
	evec2 = np.dot(evec.transpose(), SInvHalf)

	#P2 = np.zeros((nbas, nbas))
	#for i in range(0, nbas):
	#	for m in range(0, nbas):
	#		for n in range(0, nbas):
	#			P2[m, n] += evec2[i, m] * evec2[i, n] * eval[i]	

	return eval, evec2

def contravariant_basis(C, Sao):

	SC = np.dot(C, Sao)
	CtSC = np.dot(SC, C.transpose())
	metric = symmetric_matrix_power(CtSC, -1.0)
	return np.dot(metric, C)

def MO_overlap(C1, C2, Sao):
	SC2 = np.dot(C2, Sao)
	return np.dot(SC2, C1.transpose())


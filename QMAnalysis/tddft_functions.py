# -*- coding: utf-8 -*-
#
# Copyright (c), the QMAnalysis development team
#
#
#

""" CIS/TDDFT related routines  """


import sys, os
import numpy as np
from .molecule import Molecule
from .util import *
from .math_util import *

""" Co: occupied molecular orbitals, nbas-by-nocc
    Cv: occupied molecular orbitals, nbas-by-nvir
    X:  amplitudes, nvir-by-nocc """

#Rw = Cv X Co^t
def compute_transition_density_matrix(nbas, nocc, nvir, Co, Cv, X):
	#print 'nbas=', nbas, 'nocc=', nocc, 'nvir=', nvir
	Co2 = np.reshape(Co, (-1, nbas))
	Cv2 = np.reshape(Cv, (-1, nbas))
	X2  = np.reshape(X, (-1, nvir))
	CvX = np.dot(X2, Cv2)
	tdm = np.dot(Co2.transpose(), CvX).ravel()
	return tdm

#Pw = Cv X X^t Cv^t - Co X^t X Co^t = Cv X (Cv X)^t - Co X^t (Co X^t)^t 
def compute_difference_density_matrix(ddm, attdm, detdm, nbas, nocc, nvir, Co, Cv, X):

        #print('nbas=', nbas, 'nocc=', nocc, 'nvir=', nvir, len(Co), len(Cv), len(X))
        Co2  = np.reshape (Co, (-1, nbas))
        Cv2  = np.reshape(Cv,(-1, nbas))
        X2   = np.reshape(X, (-1, nvir))

        CvX  = np.dot(X2, Cv2)
        A    = np.dot(CvX.transpose(), CvX)
        attdm[:] = A.ravel()
	#matrix_print_2d2(A, 6, "Attach")

        CoXt = np.dot(X2.transpose(),Co2)
        D    = np.dot( CoXt.transpose(),CoXt)
        detdm[:] = D.ravel()
	#matrix_print_2d2(D, 6, "Detach")

        ddm[:] =  A.ravel()-D.ravel()
	#matrix_print_2(ddm, nbas, nbas, 6, "difference density matrix")

        return
        
        

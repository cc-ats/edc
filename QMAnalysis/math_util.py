import sys, string
write = sys.stdout.write
import numpy as np
from scipy import linalg
from math import sqrt, pow

def matrix_print_1d(array, m, n, ncols, title):
	""" printing a rectangular matrix, ncols columns per batch """

	write(title+':\n')
	nbatches = int(n/ncols)
	if nbatches * ncols < n: nbatches += 1
	for k in range(0, nbatches):
		write('    ')  
		j1 = ncols*k
		j2 = ncols*(k+1)
		if k == nbatches-1: j2 = n 
		for j in range(j1, j2):
			write('   %7d  ' % (j+1))
		write('\n')
		for i in range(0, m): 
			write(' %2d -' % (i+1))
			for j in range(j1, j2):
				write(' %11.7f' % array[i+j*m])
			write('\n')

def matrix_print_2d(array, ncols, title):
	""" printing a rectangular matrix, ncols columns per batch """

	write(title+'\n')
	m = array.shape[1]
	n = array.shape[0]
	#write('m=%d n=%d\n' % (m, n))
	nbatches = int(n/ncols)
	if nbatches * ncols < n: nbatches += 1
	for k in range(0, nbatches):
		write('     ')  
		j1 = ncols*k
		j2 = ncols*(k+1)
		if k == nbatches-1: j2 = n 
		for j in range(j1, j2):
			write('   %7d  ' % (j+1))
		write('\n')
		for i in range(0, m): 
			write(' %3d -' % (i+1))
			for j in range(j1, j2):
				write(' %11.6f' % array[j,i])
			write('\n')

def matrix_int_print_2d(array, m, n, title):
	write(title+":\n")
	write('   ')  
	for j in range(0, n):
		write(' %3d ' % (j+1))
	write("\n")
	for i in range(0, m):
		write('%3d' % (i+1))
		for j in range(0, n):
			write(' %3d ' % (array[i+j*m]))
		write("\n")
	write("\n")

def symmetric_matrix_sqrt(array):
	
	n = array.shape[0]
	eval, evec = linalg.eigh(array)
	print("eval: ", eval)
	evec2 = np.zeros((n,n))
	for j in range(0, n):
		if eval[j] < 0.00000001:
			scale = 0.0
		else: 
			scale = sqrt(sqrt(eval[j]))
		for i in range(0, n):	
			evec2[i,j] = evec[i,j] * scale
	prod = np.dot(evec2, evec2.transpose())
	return prod

def symmetric_matrix_power(array, alpha):
	
	n = array.shape[0]
	eval, evec = linalg.eigh(array)
	print("eval: ", eval)
	evec2 = np.zeros((n,n))
	for j in range(0, n):
		if eval[j] < 0.00000001:
			scale = 0.0
		else: 
			scale = pow(eval[j], 0.5*alpha)
		for i in range(0, n):	
			evec2[i,j] = evec[i,j] * scale
	product = np.dot(evec2, evec2.transpose())
	return product

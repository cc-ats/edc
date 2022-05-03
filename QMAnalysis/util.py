# -*- coding: utf-8 -*-
#
# Copyright (c), the QMAnalysis development team
#
#
#

""" misc utility functions """

import numpy as np
import sys, os
from math import sqrt, cos, sin

PI   = 3.141592654
BOHR = 0.529177249
IRT3 = 1/sqrt(3.)  
AU2KCAL = 627.5095

CAL2J   = 4.184

def float_str(x, ndigit=3):
	if ndigit == 3:
		return str(format(round(x, 3), ".3f"))

def is_integer(v):
    	try:     i = int(v)
    	except:  return False
    	return True

def is_float(v):
    	try:     i = float(v)
    	except:  return False
    	return True

def get_first_value(filename, phrase, pos):
	command1 = 'grep \"'+phrase+'\" ' + filename + ' | head -1 > tmp_file'
	os.system(command1)
	file = open("tmp_file", "r")
	lines = file.readlines()
	value = -1
	if len(lines) == 1:
		columns = lines[0].split()
		ncols   = len(columns)
		if pos == -1:
			value = columns[ncols-1]
		else:
			value = columns[pos]
	return value

def get_last_value(filename, phrase, pos):
	command1 = 'grep \"'+phrase+'\" ' + filename + ' | tail -1 > tmp_file'
	os.system(command1)
	file = open("tmp_file", "r")
	lines = file.readlines()
	value = -100000
	if len(lines) == 1:
		columns = lines[0].split()
		ncols   = len(columns)
		if pos == -1:
			value = columns[ncols-1]
		else:
			value = columns[pos]
	return value

def get_all_values(filename, phrase, pos):
	command1 = 'grep \"'+phrase+'\" ' + filename + ' > tmp_file'
	os.system(command1)
	file = open("tmp_file", "r")
	lines = file.readlines()
	values = []
	for line in lines:
		columns = line.split()
		ncols   = len(columns)
		if pos == -1:
			values.append(columns[ncols-1])
		else:
			values.append(columns[pos])
	return values

def get_all_values_with_index(filename, phrase, pos, index):
	command1 = 'grep -n \"'+phrase+'\" ' + filename + ' > tmp_file'
	os.system(command1)
	file = open("tmp_file", "r")
	lines = file.readlines()
	values = []
	for line in lines:
		columns = line.split()
		ncols   = len(columns)
		if pos == -1:
			values.append(columns[ncols-1])
		else:
			values.append(columns[pos])
		lineno = int(columns[0][0:columns[0].find(":")])
		index.append(lineno)
	return values

def get_first_value_2(filename, phrase1, phrase2, pos):
	""" both phrases might appear in the same line """

	command1 = 'grep \"'+phrase1+'\" ' + filename + ' | grep ' + phrase2 + ' | head -1 > tmp_file'
	os.system(command1)
	file = open("tmp_file", "r")
	lines = file.readlines()
	value = -1
	if len(lines) == 1:
		columns = lines[0].split()
		ncols   = len(columns)
		if pos == -1:
			value = columns[ncols-1]
		else:
			value = columns[pos]
	file.close()
	return value

def all_occurrences(filename, phrase):
	index = []
	command1 = 'grep -n \"'+phrase+'\" ' + filename + ' > tmp_file'
	os.system(command1)
	if os.stat("tmp_file").st_size > 0:
		file = open("tmp_file", "r")
		lines   = file.readlines()
		for k in range(0, len(lines)):
			columns = lines[k].split()
			col0 = columns[0]
			value = int(col0[0:col0.find(":")])
			index.append(value)
	return index

def find_occurrence(filename, phrase, first_or_last, return_nvalues):
	value1 = value2 = -1
	os.system('grep -n \"'+phrase+'\" ' + filename + ' > tmp_file')
	if os.stat("tmp_file").st_size > 0:
		if first_or_last == 1:
			command1 = 'grep -n \"'+phrase+'\" ' + filename + ' | head -1 > tmp_file'
		elif first_or_last == 0: 
			command1 = 'grep -n \"'+phrase+'\" ' + filename + ' | tail -1 > tmp_file'
		else: 
			command1 = 'grep -n \"'+phrase+'\" ' + filename + ' | head -'+str(first_or_last)+' | tail -1 > tmp_file'
		#print("command1: ", command1)
		os.system(command1)
		file = open("tmp_file", "r")
		lines   = file.readlines()
		columns = lines[0].split()
		col0 = columns[0]
		iend = 0
		for k in range(0, len(col0)):
			if col0[k:k+1] == ":":
				iend = k
				break
		value1 = int(col0[0:iend])
		if return_nvalues == 1: 
			value2 = int(columns[-1])
	if return_nvalues == 1:
		return [value1, value2]
	else: 
		return value1

def first_occurrence(filename, phrase, return_nvalues=0):
	return find_occurrence(filename, phrase, 1, return_nvalues)

def last_occurrence(filename, phrase, return_nvalues=0):
	return find_occurrence(filename, phrase, 0, return_nvalues)

def get_integer_block(filename, phrase):
	iline = find_occurrence(filename, phrase, 1, 0)
	outf = open(filename, "r")
	values = []
	for line in outf.readlines()[iline:iline+10]:	
		cols = line.split()
		ncol = len(cols)
		if ncol > 0 and is_integer(cols[0]):
			values.append(int(cols[0]))
		else:
			return values

# read a block of data printed using MatPrint
def get_real_data_block(filename, M, N, Words,first_occurrence=1):
	words = Words.split()
	len_words = len(words)

	array = np.zeros((N, M))
	iline = find_occurrence(filename, Words, first_occurrence, 0)
	#print("iline=", iline)

	outf = open(filename, "r")
	start = 1 
	column_indices =[0, 0, 0, 0, 0, 0, 0, 0]
	for line in outf.readlines()[iline:-1]:
		columns= line.split()
		ncols = len(columns)
		if start == 1 and ncols >= 1:
			#print('line=', line)
			if is_integer(columns[0]) == 0:    #end of block
				start = 0 
			elif ncols > 1 and is_integer(columns[1]) == 0: # real data
				row_index = int(columns[0]) - 1 
				#print('row_index=', row_index)
				for k in range(1, ncols):
					if column_indices[k-1] <= N-1:   #sometimes the dimension is wrong
						array[column_indices[k-1], row_index] = float(columns[k])
			else:
				for k in range(0, ncols):
					column_indices[k] = int(columns[k])-1
				#print('column_indices=', column_indices)
				start = 1 
	return array

def atmnum(atmsym):	
	if   atmsym == 'H'  or atmsym == 'h' or atmsym == '1':  return 1
	elif atmsym == 'He' or atmsym == 'he' or atmsym == '2': return 2
	elif atmsym == 'B'  or atmsym == 'b' or atmsym == '5':  return 5
	elif atmsym == 'C'  or atmsym == 'c' or atmsym == '6':  return 6
	elif atmsym == 'N'  or atmsym == 'n' or atmsym == '7':  return 7
	elif atmsym == 'O'  or atmsym == 'o' or atmsym == '8':  return 8
	elif atmsym == 'F'  or atmsym == 'f' or atmsym == '9':  return 9
	elif atmsym == 'Na' or atmsym == 'na' or atmsym == '11': return 11
	elif atmsym == 'S'  or atmsym == 's' or atmsym == '16':  return 16
	elif atmsym == 'Cl' or atmsym == 'cl' or atmsym == '17': return 17
	elif atmsym == 'K'  or atmsym == 'k'  or atmsym == '19':  return 19
	elif atmsym == 'Ca' or atmsym == 'ca' or atmsym == '19': return 20
	elif atmsym == 'Cu' or atmsym == 'cu' or atmsym == '29': return 29
	else: print("atom sym not supported: ", atmsym)

def atmsym(atmnum):
	if   atmnum ==  1: return "H"
	elif atmnum ==  2: return "He"
	elif atmnum ==  5: return "B"
	elif atmnum ==  6: return "C"
	elif atmnum ==  7: return "N"
	elif atmnum ==  8: return "O"
	elif atmnum ==  9: return "F"
	elif atmnum == 11: return "Na"
	elif atmnum == 16: return "S"
	elif atmnum == 17: return "Cl"
	elif atmnum == 19: return "K"
	elif atmnum == 20: return "Ca"
	elif atmnum == 29: return "Cu"
	else: return "UNK"

def tip5p_data(coords, atom1):

	#Coordinates
	o = [coords[3*atom1-3],coords[3*atom1-2],coords[3*atom1-1]] 
	h1 = [coords[3*atom1-0],coords[3*atom1+1],coords[3*atom1+2]]
	h2 = [coords[3*atom1+3],coords[3*atom1+4],coords[3*atom1+5]]
	m = 0.5 * np.add(h1,h2)
	
	#Vectors
	h1h2 = np.subtract(h2,h1)
	mo = np.subtract(o,m)
	
	#Constants
	a = 104.52   #theta,HOH
	b = 109.47   #phi,LOL
	r1 = 0.9572 / BOHR   #HO
	r2 = 0.7 / BOHR   #LO
	
	c = np.dot(mo,mo.T)
	t = (r2*cos(b/2/180*PI))/sqrt(c)   #parameter of MO equation
	n = o+t*mo   #coordinates of N
	
	B = np.cross(mo, h1h2)
	c = np.dot(B,B)
	#print("B=", B)
	#print("c=", c)
	k = (r2*sin(b/2/180*PI))/sqrt(c)   #parameter of LN equation
	
	l1 = np.add(n,k*B)   #coordinates of L1
	l2 = np.subtract(n,k*B)   #coordinates of L2

	charges = [0.241, 0.241, -0.241, -0.241]
	xyz = [h1[0], h1[1], h1[2], h2[0], h2[1], h2[2], l1[0], l1[1], l1[2], l2[0], l2[1], l2[2]]
	#print("xyz=", xyz)

	return charges, xyz


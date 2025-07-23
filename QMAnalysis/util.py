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
	if   atmsym.lower() == 'h'  or atmsym == '1':   return 1
	elif atmsym.lower() == 'he' or atmsym == '2':   return 2
	elif atmsym.lower() == 'li' or atmsym == '3':   return 3
	elif atmsym.lower() == 'be' or atmsym == '4':   return 4
	elif atmsym.lower() == 'b'  or atmsym == '5':   return 5
	elif atmsym.lower() == 'c'  or atmsym == '6':   return 6
	elif atmsym.lower() == 'n'  or atmsym == '7':   return 7
	elif atmsym.lower() == 'o'  or atmsym == '8':   return 8
	elif atmsym.lower() == 'f'  or atmsym == '9':   return 9
	elif atmsym.lower() == 'ne' or atmsym == '10':  return 10
	elif atmsym.lower() == 'na' or atmsym == '11':  return 11
	elif atmsym.lower() == 'mg' or atmsym == '12':  return 12
	elif atmsym.lower() == 'al' or atmsym == '13':  return 13
	elif atmsym.lower() == 'si' or atmsym == '14':  return 14
	elif atmsym.lower() == 'p'  or atmsym == '15':  return 15
	elif atmsym.lower() == 's'  or atmsym == '16':  return 16
	elif atmsym.lower() == 'cl' or atmsym == '17':  return 17
	elif atmsym.lower() == 'ar' or atmsym == '18':  return 18
	elif atmsym.lower() == 'k'  or atmsym == '19':  return 19
	elif atmsym.lower() == 'ca' or atmsym == '20':  return 20
	elif atmsym.lower() == 'sc' or atmsym == '21':  return 21
	elif atmsym.lower() == 'ti' or atmsym == '22':  return 22
	elif atmsym.lower() == 'v'  or atmsym == '23':  return 23
	elif atmsym.lower() == 'cr' or atmsym == '24':  return 24
	elif atmsym.lower() == 'mn' or atmsym == '25':  return 25
	elif atmsym.lower() == 'fe' or atmsym == '26':  return 26
	elif atmsym.lower() == 'co' or atmsym == '27':  return 27
	elif atmsym.lower() == 'ni' or atmsym == '28':  return 28
	elif atmsym.lower() == 'cu' or atmsym == '29':  return 29
	elif atmsym.lower() == 'zn' or atmsym == '30':  return 30
	elif atmsym.lower() == 'ga' or atmsym == '31':  return 31
	elif atmsym.lower() == 'ge' or atmsym == '32':  return 32
	elif atmsym.lower() == 'as' or atmsym == '33':  return 33
	elif atmsym.lower() == 'se' or atmsym == '34':  return 34
	elif atmsym.lower() == 'br' or atmsym == '35':  return 35
	elif atmsym.lower() == 'kr' or atmsym == '36':  return 36
	elif atmsym.lower() == 'rb' or atmsym == '37':  return 37
	elif atmsym.lower() == 'sr' or atmsym == '38':  return 38
	elif atmsym.lower() == 'y'  or atmsym == '39':  return 39
	elif atmsym.lower() == 'zr' or atmsym == '40':  return 40
	elif atmsym.lower() == 'nb' or atmsym == '41':  return 41
	elif atmsym.lower() == 'mo' or atmsym == '42':  return 42
	elif atmsym.lower() == 'tc' or atmsym == '43':  return 43
	elif atmsym.lower() == 'ru' or atmsym == '44':  return 44
	elif atmsym.lower() == 'rh' or atmsym == '45':  return 45
	elif atmsym.lower() == 'pd' or atmsym == '46':  return 46
	elif atmsym.lower() == 'ag' or atmsym == '47':  return 47
	elif atmsym.lower() == 'cd' or atmsym == '48':  return 48
	elif atmsym.lower() == 'in' or atmsym == '49':  return 49
	elif atmsym.lower() == 'sn' or atmsym == '50':  return 50
	elif atmsym.lower() == 'sb' or atmsym == '51':  return 51
	elif atmsym.lower() == 'te' or atmsym == '52':  return 52
	elif atmsym.lower() == 'i'  or atmsym == '53':  return 53
	elif atmsym.lower() == 'xe' or atmsym == '54':  return 54
	elif atmsym.lower() == 'cs' or atmsym == '55':  return 55
	elif atmsym.lower() == 'ba' or atmsym == '56':  return 56
	elif atmsym.lower() == 'la' or atmsym == '57':  return 57
	elif atmsym.lower() == 'ce' or atmsym == '58':  return 58
	elif atmsym.lower() == 'pr' or atmsym == '59':  return 59
	elif atmsym.lower() == 'nd' or atmsym == '60':  return 60
	elif atmsym.lower() == 'pm' or atmsym == '61':  return 61
	elif atmsym.lower() == 'sm' or atmsym == '62':  return 62
	elif atmsym.lower() == 'eu' or atmsym == '63':  return 63
	elif atmsym.lower() == 'gd' or atmsym == '64':  return 64
	elif atmsym.lower() == 'tb' or atmsym == '65':  return 65
	elif atmsym.lower() == 'dy' or atmsym == '66':  return 66
	elif atmsym.lower() == 'ho' or atmsym == '67':  return 67
	elif atmsym.lower() == 'er' or atmsym == '68':  return 68
	elif atmsym.lower() == 'rm' or atmsym == '69':  return 69
	elif atmsym.lower() == 'yb' or atmsym == '70':  return 70
	elif atmsym.lower() == 'lu' or atmsym == '71':  return 71
	elif atmsym.lower() == 'hf' or atmsym == '72':  return 72
	elif atmsym.lower() == 'ta' or atmsym == '73':  return 73
	elif atmsym.lower() == 'w'  or atmsym == '74':  return 74
	elif atmsym.lower() == 're' or atmsym == '75':  return 75
	elif atmsym.lower() == 'os' or atmsym == '76':  return 76
	elif atmsym.lower() == 'ir' or atmsym == '77':  return 77
	elif atmsym.lower() == 'pt' or atmsym == '78':  return 78
	elif atmsym.lower() == 'au' or atmsym == '79':  return 79
	elif atmsym.lower() == 'hg' or atmsym == '80':  return 80
	elif atmsym.lower() == 'tl' or atmsym == '81':  return 81
	elif atmsym.lower() == 'pb' or atmsym == '82':  return 82
	elif atmsym.lower() == 'bi' or atmsym == '83':  return 83
	elif atmsym.lower() == 'po' or atmsym == '84':  return 84
	elif atmsym.lower() == 'at' or atmsym == '85':  return 85
	elif atmsym.lower() == 'rn' or atmsym == '86':  return 86
	elif atmsym.lower() == 'fr' or atmsym == '87':  return 87
	elif atmsym.lower() == 'ra' or atmsym == '88':  return 88
	elif atmsym.lower() == 'ac' or atmsym == '89':  return 89
	elif atmsym.lower() == 'th' or atmsym == '90':  return 90
	elif atmsym.lower() == 'pa' or atmsym == '91':  return 91
	elif atmsym.lower() == 'u'  or atmsym == '92':  return 92
	elif atmsym.lower() == 'np' or atmsym == '93':  return 93
	elif atmsym.lower() == 'pu' or atmsym == '94':  return 94
	elif atmsym.lower() == 'am' or atmsym == '95':  return 95
	elif atmsym.lower() == 'cm' or atmsym == '96':  return 96
	elif atmsym.lower() == 'bk' or atmsym == '97':  return 97
	elif atmsym.lower() == 'cf' or atmsym == '98':  return 98
	elif atmsym.lower() == 'es' or atmsym == '99':  return 99
	elif atmsym.lower() == 'fm' or atmsym == '100':  return 100
	elif atmsym.lower() == 'md' or atmsym == '101':  return 101
	elif atmsym.lower() == 'no' or atmsym == '102':  return 102
	elif atmsym.lower() == 'lr' or atmsym == '103':  return 103
	elif atmsym.lower() == 'rf' or atmsym == '104':  return 104
	elif atmsym.lower() == 'db' or atmsym == '105':  return 105
	elif atmsym.lower() == 'sg' or atmsym == '106':  return 106
	elif atmsym.lower() == 'bh' or atmsym == '107':  return 107
	elif atmsym.lower() == 'hs' or atmsym == '108':  return 108
	elif atmsym.lower() == 'mt' or atmsym == '109':  return 109
	elif atmsym.lower() == 'ds' or atmsym == '110':  return 110
	elif atmsym.lower() == 'rg' or atmsym == '111':  return 111
	elif atmsym.lower() == 'cn' or atmsym == '112':  return 112
	elif atmsym.lower() == 'nh' or atmsym == '113':  return 113
	elif atmsym.lower() == 'fl' or atmsym == '114':  return 114
	elif atmsym.lower() == 'mc' or atmsym == '115':  return 115
	elif atmsym.lower() == 'lv' or atmsym == '116':  return 116
	elif atmsym.lower() == 'ts' or atmsym == '117':  return 117
	elif atmsym.lower() == 'og' or atmsym == '118':  return 118
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

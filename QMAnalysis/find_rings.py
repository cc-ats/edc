import string, sys
from math import sqrt

def find_5_membered_rings(connectivity):
	natoms = int(sqrt(len(connectivity)))
	print('natoms=', natoms)
	ring_5_index = []
	for k in range(0, natoms):
		ring_5_index.append(0)
	count = 0
	for i in range(0, natoms):
		for j in range(0, natoms):
			if connectivity[j+i*natoms] == 2:
				for k in range(0, natoms):
					if k != j and connectivity[k+i*natoms] == 2 and connectivity[j+k*natoms] == 1:
						#print i+1, j+1, k+1
						if ring_5_index[i] + ring_5_index[j] + ring_5_index[k] == 0:
							count += 1
							ring_5_index[i] = ring_5_index[j] = ring_5_index[k] = count
						elif ring_5_index[i] != 0:
							if ring_5_index[j] == 0 : ring_5_index[j] = ring_5_index[i]
							if ring_5_index[k] == 0 : ring_5_index[k] = ring_5_index[i]
						elif ring_5_index[j] != 0:
							if ring_5_index[i] == 0 : ring_5_index[i] = ring_5_index[j]
							if ring_5_index[k] == 0 : ring_5_index[k] = ring_5_index[j]
						elif ring_5_index[k] != 0:
							if ring_5_index[i] == 0 : ring_5_index[i] = ring_5_index[k]
							if ring_5_index[j] == 0 : ring_5_index[j] = ring_5_index[k]

	return ring_5_index
	#for k in range(0, natoms): print k+1, ring_5_index[k]
							

def find_6_membered_rings(connectivity):
	natoms = int(sqrt(len(connectivity)))
	#print 'natoms=', natoms
	ring_6_index = []
	for k in range(0, natoms):
		ring_6_index.append(0)
	count = 0
	for i in range(0, natoms):
		for j in range(0, natoms):
			if connectivity[j+i*natoms] == 3:
				count2 = 0
				index = []
				for k in range(0, natoms):
					if k != i and k != j and connectivity[k+j*natoms] == 1 and connectivity[k+i*natoms] == 2:
						index.append(k)
						count2 += 1
				if count2 == 2: 
					k = index[0]
					l = index[1]
					if ring_6_index[i] + ring_6_index[j] + ring_6_index[k] + ring_6_index[l] == 0:
						count += 1
						ring_6_index[i] = ring_6_index[j] = ring_6_index[k] = ring_6_index[l] = count;
					elif ring_6_index[i] != 0:
						if ring_6_index[j] == 0 : ring_6_index[j] = ring_6_index[i]
						if ring_6_index[k] == 0 : ring_6_index[k] = ring_6_index[i]
						if ring_6_index[l] == 0 : ring_6_index[l] = ring_6_index[i]
					elif ring_6_index[j] != 0:
						if ring_6_index[i] == 0 : ring_6_index[i] = ring_6_index[j]
						if ring_6_index[k] == 0 : ring_6_index[k] = ring_6_index[j]
						if ring_6_index[l] == 0 : ring_6_index[l] = ring_6_index[j]
					elif ring_6_index[k] != 0:
						if ring_6_index[i] == 0 : ring_6_index[i] = ring_6_index[k]
						if ring_6_index[j] == 0 : ring_6_index[j] = ring_6_index[k]
						if ring_6_index[l] == 0 : ring_6_index[l] = ring_6_index[k]
					elif ring_6_index[l] != 0:
						if ring_6_index[i] == 0 : ring_6_index[i] = ring_6_index[l]
						if ring_6_index[j] == 0 : ring_6_index[j] = ring_6_index[l]
						if ring_6_index[k] == 0 : ring_6_index[k] = ring_6_index[l]

	#for k in range(0, natoms): print k+1, ring_6_index[k]
	return ring_6_index



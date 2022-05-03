import numpy as np
import string, sys
import sys
write = sys.stdout.write
from .find_rings import *

def find_thiazole_atoms(atomsym, atomnum, connectivity):

	natoms = len(atomsym)
	ring_5_index = find_5_membered_rings(connectivity)
	ring_6_index = find_6_membered_rings(connectivity)
	print('ring_5_index=', ring_5_index)
	print('ring_6_index=', ring_6_index)

	n_5_rings = 0
	for k in range(0, natoms):
		if ring_5_index[k] > n_5_rings:
			n_5_rings = ring_5_index[k]

	print('n_5_rings=', n_5_rings)
	is_thiazole_ring = []
	for k in range(0, n_5_rings):
		is_thiazole_ring.append(1)

	for k in range(0, natoms):
		if ring_5_index[k] > 0 and ring_6_index[k] > 0:
			is_thiazole_ring[ring_5_index[k]-1] = 0
	print('is_thiazole_ring=', is_thiazole_ring)

	#find the bridge carbon
	atom_carbon    = -1
	atom_neighbor  = -1
	for k in range(0, natoms):
		if atomsym[k] == 'C' and ring_5_index[k] > 0 and is_thiazole_ring[ring_5_index[k]-1] == 1:
			#print 'find carbon', k+1
			nitrogen_neighbor = -1
			sulfur_neighbor = -1
			for m in range(0, natoms):
				if m != k and connectivity[m+k*natoms] == 1 and ring_5_index[m] == ring_5_index[k]:
					#print m, atomsym[m]
					if atomsym[m] == 'N':
						nitrogen_neighbor = m
					elif atomsym[m] == 'O' or atomsym[m] == 'S' or atomsym[m] == 'Se':
						sulfur_neighbor = m
			if nitrogen_neighbor != -1 and sulfur_neighbor != -1:
					#print k+1, nitrogen_neighbor+1, sulfur_neighbor+1
					for m in range(0, natoms):
						if m!= k and connectivity[m+k*natoms] == 1 and ring_5_index[m] != ring_5_index[k]:
							#print 'find neighbor', m+1
							atom_carbon = k
							atom_neighbor = m

	print('atom_carbon=', atom_carbon+1)
	print('atom_neighbor=', atom_neighbor+1)

	on_thiazole_ring = []
	for k in range(0, natoms):
		if ring_5_index[k] > 0 and is_thiazole_ring[ring_5_index[k]-1] == 1:
			on_thiazole_ring.append(1)
		else:
			on_thiazole_ring.append(0)

	for k in range(0, 5):
		for i in range(0, natoms):
			if on_thiazole_ring[i] == 1:
				for j in range(0, natoms):
					if  j != i and j != atom_neighbor and connectivity[j+i*natoms] == 1:
						on_thiazole_ring[j] = 1

	write("thiazole atoms: ")
	for k in range(0, natoms):
		if on_thiazole_ring[k] > 0:
			write("%3d " % (k+1))
	write("\n")

	return on_thiazole_ring


from .molecule import Molecule
from .qm_system import QMSystem
from .cubic_grid import CubicGrid
from .xyz_file_parser import read_all_geometries, compare_molecules
from .util import is_integer, atmsym
import sys, time
write = sys.stdout.write
import numpy as np
from scipy import interpolate

def edc_set_up_structures(fchk1, fchk2, ircfile):

    qm1 = QMSystem()
    qm1.setup_from_qchem_fchk(fchk1)
    mol1 = qm1.molecule

    qm2 = QMSystem()
    qm2.setup_from_qchem_fchk(fchk2)
    mol2 = qm2.molecule

    frames = read_all_geometries(ircfile, qm1.natoms)

    if compare_molecules(frames[0], mol1) and compare_molecules(frames[-1], mol2):
        forward = 1
        print("IRC file is found to go from mol1 to mol2")
    elif compare_molecules(frames[-1], mol1) and compare_molecules(frames[0], mol2):
        forward = 0
        print("IRC file is found to go from mol2 to mol1")
    else:
        print("Wrong IRC file")
        exit(-1)

    return qm1, qm2, frames, forward

def edc_atomic_density_fits(mol):

    elements = []
    density_fits = []
    for element in mol.atmsym:
        if not (element in elements):
            elements.append(element)
            radfile = "/home/yihan/GitHub/QMAnalysis/src/qm_analysis/data/"+element+".Rad"
            rad, den = np.loadtxt(radfile, unpack='true')
            fit =  interpolate.splrep(rad, den, s=0)
            density_fits.append(fit)
    return elements, density_fits


def compute_atomic_densities(dists, mol, elements_fits, fits, ngrid):
    natoms = mol.natoms
    values = np.zeros((natoms,ngrid))
    for iatom in range(0, natoms):
        iElement = elements_fits.index(mol.atmsym[iatom])
        values[iatom,:] = interpolate.splev(dists[iatom,:], fits[iElement])
    return values

def grid_extrapolation(frames, forward, elements_fit, atomic_density_fits, gridsize):

    natoms = frames[0].natoms
    nframes = len(frames)
    if forward == 1: frame1, frame2, inc = 0, nframes-1, 1
    else: frame1, frame2, inc = nframes-1, 0, -1

    #initial grid
    grid = CubicGrid()
    grid.setup(frames[frame1], 2.0, gridsize)
    ngrid = grid.ngrid
    grid_new = np.reshape(grid.pos, (ngrid, 3)).transpose()
    print("grid_new shape", grid_new.shape)

    for m in range(frame1, frame2, inc):
        write("Extrapolating to the %3d-th frame, " %(m))
        timem = time.time()
        dx = np.zeros(ngrid)
        dy = np.zeros(ngrid)
        dz = np.zeros(ngrid)
        grid_atom_dist = np.zeros((natoms, ngrid))
        for iAtom in range(0, natoms):
            dx[:] = frames[m].coords[iAtom*3 + 0] - grid_new[0,:]
            dy[:] = frames[m].coords[iAtom*3 + 1] - grid_new[1,:]
            dz[:] = frames[m].coords[iAtom*3 + 2] - grid_new[2,:]
            grid_atom_dist[iAtom, :] = np.sqrt(dx*dx + dy*dy + dz*dz)
        atomic_density = compute_atomic_densities(grid_atom_dist, frames[0], elements_fit, atomic_density_fits, ngrid)
        total_density = np.sum(atomic_density, axis=0)
        weights = np.zeros((natoms, ngrid))
        for iAtom in range(0, natoms):
            weights[iAtom,:] = atomic_density[iAtom,:] / total_density[:]
        for iAtom in range(0, natoms):
            for i in range(0, 3):
                grid_new[i,:] += weights[iAtom,:] * (frames[m+inc].coords[3*iAtom+i] - frames[m].coords[3*iAtom+i])
        write("  Time: %.2fs\n" %(round(- timem + time.time(),2)))

    return grid_new.transpose().ravel()

def plot_density_difference(qm1, qm2, forward, distorted_grid, gridsize, plot_name):

    qm1.plot_total_density(gridsize)
    qm2.user_defined_grid(distorted_grid)
    qm2.plot_total_density()
    difference_density_on_grid = np.subtract(qm2.density_on_grid, qm1.density_on_grid)
    qm1.plot_grid_data(plot_name, difference_density_on_grid)


def edc_analysis(fchk1, fchk2, ircfile, gridsize, plot_name='unk'):

    qm1, qm2, frames, forward = edc_set_up_structures(fchk1, fchk2, ircfile)
    elements_fits, atomic_density_fits = edc_atomic_density_fits(frames[0])
    distorted_grid = grid_extrapolation(frames, forward, elements_fits, atomic_density_fits, gridsize)
    plot_density_difference(qm1, qm2, forward, distorted_grid, gridsize, plot_name)



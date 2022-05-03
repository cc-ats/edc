import sys, os, time, glob
import numpy as np
import QMAnalysis.src.qm_analysis as qm
write = sys.stdout.write

if len(sys.argv) > 1:
    direction = sys.argv[1]
else:
    write("Provide a direction for the movie:\n")
    write("1 for forward\n")
    write("-1 for reverse\n")
    exit(-1)
gridsize = 0.1
xyzfile_name = 'tmp_irc.xyz'
filelist = glob.glob('*.fchk')
if int(direction) == 1:
    filelist.sort()
    movie_dir = 'f'
elif int(direction) == -1:
    filelist.sort(reverse=True)
    movie_dir = 'f'

for index, file in enumerate(filelist):
    if index+1 < len(filelist):
        write("Calculating frame %d to frame %d" %(index, index+1))
        fchk_i = filelist[index]
        fchk_f = filelist[index+1]
        print(fchk_i)
        print(fchk_f)
        fchk1 = qm.QChemFChk(fchk_i)
        fchk2 = qm.QChemFChk(fchk_f)
        geom1 = fchk1.geometry
        geom2 = fchk2.geometry
        if index < 10:
            jobname = movie_dir+"0"+str(index)+"_frame"
        else:
            jobname = movie_dir+str(index)+"_frame"
        with open(xyzfile_name, 'w') as xyzfile:
            xyzfile.write(str(geom1.natoms)+"\n")
            xyzfile.write(" \n")
            for atom in range(geom1.natoms):
                xyzfile.write("{:<3}".format(geom1.atmsym[atom]))
                xyzfile.write("{:>12.8f}".format(geom1.coords[3*atom]*0.529177249))
                xyzfile.write("{:>12.8f}".format(geom1.coords[3*atom+1]*0.529177249))
                xyzfile.write("{:>12.8f}".format(geom1.coords[3*atom+2]*0.529177249)+"\n")
            xyzfile.write(str(geom2.natoms)+"\n")
            xyzfile.write(" \n")
            for atom in range(geom2.natoms):
                xyzfile.write("{:<3}".format(geom2.atmsym[atom]))
                xyzfile.write("{:>12.8f}".format(geom2.coords[3*atom]*0.529177249))
                xyzfile.write("{:>12.8f}".format(geom2.coords[3*atom+1]*0.529177249))
                xyzfile.write("{:>12.8f}".format(geom2.coords[3*atom+2]*0.529177249)+"\n")
        qm.edc_analysis(fchk_i, fchk_f, xyzfile_name, gridsize, jobname)

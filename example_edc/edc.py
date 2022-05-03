import sys 
import QMAnalysis.src.qm_analysis as qm
write = sys.stdout.write

fchk1 = sys.argv[1]
fchk2 = sys.argv[2]
ircfile = sys.argv[3]
if len(sys.argv) > 4: 
    plot_name = sys.argv[4]
else: 
    print("Give a plot name")
    exit(-1)

gridsize = 0.1
qm.edc_analysis(fchk1, fchk2, ircfile, gridsize, plot_name)

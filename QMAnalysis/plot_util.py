import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def scatter_plot(v1, v2, xlabel, ylabel, figfile):

	#print 'length=', len(v1)
	#print 'v1=', v1
	#print 'v2=', v2
	plt.plot(v1, v2, marker='.', markersize=8, linestyle="None")
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.savefig(figfile)
	plt.close()

def add_arrow(fig, x0, y0, dx, dy, color, Fill):
	ax = fig.add_subplot(111)
	ax.add_patch(
	patches.Arrow(
		x0, y0, dx, dy, 
		width=1,      # width
		color=color,
		fill=Fill
	)   
)


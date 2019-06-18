#!/usr/bin/env python

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys

#mpl.rcParams['pdf.fonttype'] = 42

#colors = ['red', 'blue', 'green', 'fuchsia', 'darkorange', 'aqua', 'dodgerblue', 'goldenrod', 'lightseagreen', 'darkviolet']
colors = ['red', 'blue', 'green', 'darkviolet', 'lightseagreen', 'orange']
markers = ['s', 'o', '^', '+', 'x', 'd']

multList = ['4', '5', '6', '10', '100', '1000']
nLoopsList = ['100', '10000', '1000000']
binwidthList = ['0.5', '0.1', '0.05', '0.025', '0.01', '0.005']
#binwidthList = ['0.05', '0.025', '0.01', '0.005']

#=====================================
def compare_binwidths(mult, nLoops):
	plotfontsize = 16
	f, ax = plt.subplots(1)
	
	lw = 2

	xpts = np.linspace(-0.2, 0.2, 1001.0)
	exactCurve1 = 1+np.exp((-4.5/0.19733**2+1.0/(2.0*0.5**2))*xpts**2)
	plt.plot(xpts, exactCurve1, color='black', linewidth=1, label='Exact result')

	color_idx = dict(zip(binwidthList, (np.linspace(0.0, 1.0, len(binwidthList))).tolist()))
	
	count = 0
	for bw in binwidthList:
		fileToPlot = "./auto/results_mult" + str(mult) + "_nLoops" + str(nLoops) + "_bw" + bw + "/HBT_pipiCF.dat"
		data = np.loadtxt(fileToPlot)
		# in center cell, just take error to be zero
		data[np.where(np.abs(data[:,5])<1.0e-10),10] = 0.0
		#thisColor = plt.cm.hsv(color_idx[bw])
		thisColor = colors[count]
		thisMarker = markers[count]
		count += 1
		
		plt.plot(data[:,5], data[:,9], color=thisColor, markersize=0)
		plt.errorbar(data[:,5], data[:,9], yerr=data[:,10], color=thisColor, markersize=0, label='Event generator', fmt='o')
		
		# include error bands
		plt.fill_between(data[:,5],\
							data[:,9]-data[:,10],\
							data[:,9]+data[:,10],\
							alpha=0.05, edgecolor=thisColor,\
							facecolor=thisColor)


	
	# plot a unity line to guide the eye
	plt.plot(xpts, 1.0 + 0.0*xpts, color='black', linewidth=1)
	
	# label stuff
	plt.xlabel(r'$Q$ (GeV)', fontsize=20)
	plt.ylabel(r'$C(Q, \vec{K} = \vec{0})$', fontsize=20)
	plt.xlim( -0.2, 0.2 )
	plt.ylim( 0.75, 2.25 )
	#plt.legend(loc='upper center', fontsize=20)

	# last CMD line argument is output filename
	#plt.savefig(sys.argv[2])
	plt.show()


def exactCurve(x):
	return 1+np.exp((-4.5/0.19733**2+1.0/(2.0*0.5**2))*x**2)

#=====================================
def relative_compare_binwidths(mult, nLoops):
	plotfontsize = 16
	f, ax = plt.subplots(1)
	
	lw = 2

	xpts = np.linspace(-0.2, 0.2, 1001.0)
	#exactCurve1 = 1+np.exp((-4.5/0.19733**2+1.0/(2.0*0.5**2))*xpts**2)
	#plt.plot(xpts, exactCurve1, color='black', linewidth=1, label='Exact result')

	color_idx = dict(zip(binwidthList,\
						( np.arange( 0.0, 1.0, 1.0/( 1.0 + len( binwidthList ) ) )[1:] ).tolist()\
					) )

	count = 0
	for bw in binwidthList:
		fileToPlot = "./auto/results_mult" + str(mult) + "_nLoops" + str(nLoops) + "_bw" + bw + "/HBT_pipiCF.dat"
		data = np.loadtxt(fileToPlot)
		# in center cell, just take error to be zero
		data[np.where(np.abs(data[:,5])<1.0e-10),10] = 0.0
		exactPoints = exactCurve(data[:,5])
		
		#thisColor = plt.cm.hsv(color_idx[bw])
		thisColor = colors[count]
		thisMarker = markers[count]
		count += 1
		
		#plt.plot(data[:,5], data[:,9]-exactPoints, color=thisColor, markersize=10, marker=thisMarker)
		#plt.errorbar(data[:,5], data[:,9]-exactPoints, yerr=data[:,10], color=thisColor, markersize=0, label='Event generator')
		plt.plot(data[:,5], (data[:,9]-exactPoints)/(data[:,10]+1.e-10), color=thisColor, markersize=10, marker=thisMarker)
		plt.errorbar(data[:,5], (data[:,9]-exactPoints)/(data[:,10]+1.e-10), yerr=1.0, color=thisColor, markersize=0, label='Event generator')
				
		# include error bands
		'''plt.fill_between(data[:,5],\
							data[:,9]-exactPoints-data[:,10],\
							data[:,9]-exactPoints+data[:,10],\
							alpha=0.05, edgecolor=thisColor,\
							facecolor=thisColor)'''
		plt.fill_between(data[:,5],\
							(data[:,9]-exactPoints-data[:,10])/(data[:,10]+1.e-10),\
							(data[:,9]-exactPoints+data[:,10])/(data[:,10]+1.e-10),\
							alpha=0.05, edgecolor=thisColor,\
							facecolor=thisColor)
	
	# plot a unity line to guide the eye
	plt.plot(xpts, 0.0*xpts, color='black', linewidth=1)
	
	# label stuff
	plt.xlabel(r'$Q$ (GeV)', fontsize=20)
	plt.ylabel(r'$C(Q, \vec{K} = \vec{0})$', fontsize=20)
	plt.xlim( -0.2, 0.2 )
	plt.ylim( -5.0, 5.0 )
	#plt.legend(loc='upper center', fontsize=20)

	# last CMD line argument is output filename
	#plt.savefig(sys.argv[2])
	plt.show()

#=====================================

if __name__ == "__main__":
	#compare_binwidths(5, 100)
	relative_compare_binwidths(5, 100)
	print 'Finished all.'

#End of file

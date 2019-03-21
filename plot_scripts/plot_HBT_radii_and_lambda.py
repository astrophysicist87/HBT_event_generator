import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

# list of files
files = sys.argv[1:]

# plots styles, etc.
colors = ['blue', 'red', 'green', 'purple', 'orange', 'yellow']

centralities = ['0-10%', '10-40%', '40-100%']

colsToPlot = { 'HBT': 4 , 'lambda': 3 }
yLabelsToPlot = { 'HBT': r'$R^2 (\mathrm{fm}^2)$' , 'lambda': r'$\lambda$' }

def make_plot(objToPlot):
	plotfontsize = 16
	f, ax = plt.subplots(1)

	theFigureTitle = 'pp at 13 TeV ($10^7$ events)'
	f.suptitle(theFigureTitle, fontsize=plotfontsize+5)

	col = colsToPlot[objToPlot]
	yLabel = yLabelsToPlot[objToPlot]

	for iFile in xrange(len(files)):
		thisFile = files[iFile]
		thisColor = colors[iFile]
		thisCentrality = centralities[iFile]
		data = np.loadtxt(thisFile)
		data = data[( abs(data[:,1])<1.e-6 ) & ( data[:,2]<1.e-6 )]
		
		ax.plot(data[:,0], 0.0*data[:,col], color='black', linestyle='-')
		ax.plot(data[:,0], data[:,col], linestyle='-', linewidth=2, color=thisColor, label=thisCentrality)
		
		# include error bands
		plt.fill_between(data[:,0],\
							data[:,col]-data[:,col+2],\
							data[:,col]+data[:,col+2],\
							alpha=0.25, edgecolor=thisColor,\
							facecolor=thisColor)
    
	lims = [0.0, 1.0, 0.0, 1.5]
	plt.axis(lims)
	ax.set_xlabel(r'$K_T$ (GeV)', fontsize=plotfontsize)
	ax.set_ylabel(yLabel, fontsize=plotfontsize)
	ax.legend(loc='best', fontsize=plotfontsize+5)

	#miscAxesText = r'$K_T = 100$ MeV, $\Phi_K = 0$, $K_L = 0$'
	#plt.text(0.175, 1.075, miscAxesText, fontsize=plotfontsize+3)
    
	plt.show()
	outfilename = './pp_at_13TeV_KT100MeV_centrality_comparison.pdf'
	#plt.savefig(outfilename)
	print 'Saving to', outfilename
    
	plt.close()
    
if __name__ == "__main__":
    #make_plot('HBT')
    make_plot('lambda')
    print 'Finished all.'

#End of file

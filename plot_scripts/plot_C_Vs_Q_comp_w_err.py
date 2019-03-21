import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

# list of files
files = sys.argv[1:]

# plots styles, etc.
colors = ['blue', 'red', 'green', 'purple', 'orange', 'yellow']

centralities = ['0-10%', '10-40%', '40-100%']

def make_plot():
	plotfontsize = 16
	f, ax = plt.subplots(1)

	theFigureTitle = 'pp at 13 TeV ($10^7$ events)'
	f.suptitle(theFigureTitle, fontsize=plotfontsize+5)

	for iFile in xrange(len(files)):
		thisFile = files[iFile]
		thisColor = colors[iFile]
		thisCentrality = centralities[iFile]
		data = np.loadtxt(thisFile)
		data = data[( abs(data[:,1])<1.e-6 ) & ( data[:,0]==0.1 ) & ( data[:,3]>1.e-6 )]
		
		# note: correct for this in future versions of HBT code
		ax.plot(data[:,3], 0.0*data[:,7]+1, color='black', linestyle='-')
		ax.plot(data[:,3], data[:,7], linestyle='-', linewidth=2, color=thisColor, label=thisCentrality)
		
		# include error bands
		plt.fill_between(data[:,3],\
							data[:,7]-data[:,8],\
							data[:,7]+data[:,8],\
							alpha=0.25, edgecolor=thisColor,\
							facecolor=thisColor)
    
	lims = [0.05, 0.3, 0.99, 1.2]
	plt.axis(lims)
	ax.set_xlabel(r'$Q$ (GeV)', fontsize=plotfontsize)
	ax.set_ylabel(r'$C(Q)$', fontsize=plotfontsize)
	ax.legend(loc='best', fontsize=plotfontsize+5)

	miscAxesText = r'$K_T = 100$ MeV, $\Phi_K = 0$, $K_L = 0$'
	plt.text(0.175, 1.075, miscAxesText, fontsize=plotfontsize+3)
    
	plt.show()
	outfilename = './pp_at_13TeV_KT100MeV_centrality_comparison.pdf'
	#plt.savefig(outfilename)
	print 'Saving to', outfilename
    
	plt.close()
    
if __name__ == "__main__":
    make_plot()
    print 'Finished all.'

#End of file

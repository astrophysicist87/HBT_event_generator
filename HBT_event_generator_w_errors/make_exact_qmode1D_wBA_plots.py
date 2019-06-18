#!/usr/bin/env python

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import sys

mpl.rcParams['pdf.fonttype'] = 42

# files to plot come after command name and do not include last argument
fileToPlot = sys.argv[1]
exactFileToPlot = sys.argv[2]

lw = 2

xpts = np.linspace(-0.5,0.5,11)
data = np.loadtxt(fileToPlot)
exactData = np.loadtxt(exactFileToPlot)

#plt.errorbar(data[:,3], data[:,7], yerr=data[:,8], color='orange', markersize=3, fmt='o')
plt.scatter(data[:,3], data[:,4]/np.abs(data[:,3]), color='purple', label=r'$\epsilon = 1$ MeV')
plt.plot(exactData[:,0], exactData[:,1], color='purple', linewidth=1, label='Exact result')

plt.plot(xpts, 0.0*xpts, color='black', linewidth=1)

plt.xlabel(r'$Q$ (GeV)', fontsize=20)
plt.ylabel(r'$\mathrm{Num}_{\mathrm{avg}}(Q)$ (normalized)', fontsize=20)
plt.xlim( -0.5, 0.5 )
plt.ylim( 0.675, 1.025 )
plt.legend(loc=0, fontsize=20)

# last CMD line argument is output filename
plt.savefig(sys.argv[3])
#plt.show()

#End of file

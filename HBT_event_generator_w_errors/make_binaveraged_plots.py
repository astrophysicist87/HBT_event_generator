#!/usr/bin/env python

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import sys

mpl.rcParams['pdf.fonttype'] = 42

colors = ['red', 'blue', 'green']
labels = [r'Binwidth = 5 MeV', r'Binwidth = 100 MeV', r'Binwidth = 500 MeV']
normalizations = [21.0/1668568.0, 21.0/33371360.0, 21.0/166856800.0]
bValues = [1.0/0.005, 1.0/0.1, 1.0/0.5]

# files to plot come after command name and do not include last argument
filesToPlot = sys.argv[1:-1]

lw = 2

xpts = np.linspace(-0.5, 0.5, 1001.0)
#exactCurve1 = np.exp(-50.0*xpts**2)
#plt.plot(xpts, exactCurve1, color='purple', linewidth=1)

#calculated for R = 10 fm, alpha = 0.02 GeV/fm, and b = 2 1/GeV ( == 1/(0.5 GeV) )
def exactCurve2Func(b, x):
	return 0.1 * b * np.sqrt(0.5*np.pi) \
				* ( sp.erf( 5.0*(1.0-2.0*b*np.abs(x))/(b*np.sqrt(2.0)) ) \
					+ sp.erf( 5.0*(1.0+2.0*b*np.abs(x))/(b*np.sqrt(2.0)) ) )

i = 0
for filename in filesToPlot:
	data = np.loadtxt(filename)
	#plt.plot(data[:,2], data[:,6], color=colors[i], linewidth=lw, label=labels[i])
	plt.scatter(data[:,2], data[:,6]*normalizations[i], color=colors[i], label=labels[i])
	plt.plot(xpts, exactCurve2Func(bValues[i], xpts), color=colors[i], linewidth=1)
	i+=1


plt.plot(xpts, 0.0*xpts, color='black', linewidth=1)

plt.xlabel(r'$K_z$ (GeV)', fontsize=20)
plt.ylabel(r'$\mathrm{Num}(q_z, K_z)$ (normalized)', fontsize=20)
plt.xlim( -0.51, 0.51 )
plt.ylim( -0.01, 1.01 )
plt.legend(loc=0, fontsize=12)

# last CMD line argument is output filename
plt.savefig(sys.argv[-1])
#plt.show()

#End of file

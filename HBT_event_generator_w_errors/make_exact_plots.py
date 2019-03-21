#!/usr/bin/env python

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import sys

mpl.rcParams['pdf.fonttype'] = 42

colors = ['red', 'blue', 'green']

# files to plot come after command name and do not include last argument
fileToPlot = sys.argv[1]

lw = 2

xpts = np.linspace(-0.4, 0.4, 1001.0)
exactCurve1 = 1.0+np.exp(12.5*xpts**2)
plt.plot(xpts, exactCurve1, color='purple', linewidth=1, label='Exact result')


data = np.loadtxt(fileToPlot)
#plt.scatter(data[:,5], data[:,9]*normalizations[i], color=colors[i], label=labels[i])
#plt.errorbar(data[:,5], data[:,9], yerr=data[:,10], color='green', markersize=3, label='Event generator', fmt='o')
plt.errorbar(data[:,3], data[:,7], yerr=data[:,8], color='purple', markersize=3, label='Event generator', fmt='o')


plt.plot(xpts, 0.0*xpts, color='black', linewidth=1)

plt.xlabel(r'$Q$ (GeV)', fontsize=20)
plt.ylabel(r'$C(Q, \vec{K} = \vec{0})$', fontsize=20)
plt.xlim( -0.325, 0.325 )
plt.ylim( 1.8, 4.5 )
plt.legend(loc='upper center', fontsize=20)

# last CMD line argument is output filename
plt.savefig(sys.argv[2])
#plt.show()

#End of file

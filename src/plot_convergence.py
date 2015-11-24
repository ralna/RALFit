#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

try:
    problem = sys.argv[1]
except: 
    print "Please pass the problem name in by calling"
    print "   ./plot_convergence <PROBLEM_NAME> "

datatype = np.dtype({'names' : ('grad','res'),
                     'formats' : ('float','float')})

try:
    filename = "data/" + problem + "_progress.out"
    progress1 = np.loadtxt(filename+"_m1", dtype = datatype)
    progress2 = np.loadtxt(filename+"_m2", dtype = datatype)
    progress7 = np.loadtxt(filename+"_m7", dtype = datatype)
except:
    print "Sorry, no data linked with problem ", problem, " found"

plt.semilogy(progress1['grad'],label="Gauss-Newton")
plt.semilogy(progress2['grad'],label="Newton")
plt.semilogy(progress7['grad'],label="Hybrid")

plt.legend()
plt.show()

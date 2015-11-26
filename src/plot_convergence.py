#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import subprocess

short_hash = subprocess.check_output(['git','rev-parse','--short','HEAD']).strip()

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
    progress8 = np.loadtxt(filename+"_m8", dtype = datatype)
except:
    print "Sorry, no data linked with problem ", problem, " found"


all_min = np.array([ progress1['res'].min(),
                     progress2['res'].min(),
                     progress7['res'].min(),
                     progress8['res'].min()])
minvalue = all_min.min()
print 'minvalue = ', minvalue

# Plot the gradients
plt.figure(1)
plt.semilogy(progress1['grad'],label="Gauss-Newton")
plt.semilogy(progress2['grad'],label="Newton")
plt.semilogy(progress7['grad'],label="Hybrid")
plt.semilogy(progress8['grad'],label="Hybrid II")

plt.legend()
plt.title(problem+': gradients')
plt.xlabel('Iteration number')
plt.ylabel('$||J^Tr||_2$')
plt.savefig('../doc/img/'+problem+'_'+short_hash+'.png')



# Plot the residuals
plt.figure(2)
plt.semilogy(progress1['res']-minvalue,label="Gauss-Newton")
plt.semilogy(progress2['res']-minvalue,label="Newton")
plt.semilogy(progress7['res']-minvalue,label="Hybrid")
plt.semilogy(progress8['res']-minvalue,label="Hybrid II")

plt.legend()
plt.title(problem+': residuals \n minimizer = '+str(minvalue) )
plt.xlabel('Iteration number')
plt.ylabel('$1/2||r_k||^2_2 - 1/2||r_*||^2_2$')
plt.savefig('../doc/img/'+problem+'_res_'+short_hash+'.png')

plt.show()



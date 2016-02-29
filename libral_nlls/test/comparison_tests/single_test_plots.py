#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os

# Let's get the files containing the problem and control parameters from the calling command..
no_tests = len(sys.argv)-1
problem = sys.argv[1]

print "Solving problem " + problem

control_files = [None for i in range(no_tests-1)]
for i in range(no_tests-1):
    control_files[i] = sys.argv[i+2]

# read the cutest directory from the environment variables
try:
    cutestdir = os.environ.get("CUTEST")
except:
    print "ERROR: the CUTEST environment variable doesn't appear to be set"

# get the current git hash
short_hash = subprocess.check_output(['git','rev-parse','--short','HEAD']).strip()

# setup the datatype that we'll store the results in
datatype = np.dtype({'names' : ('res','grad'),
                     'formats' : ('float','float')})

# setup a python list in which to store the arrays...
progress = [None for i in range(no_tests-1)]
# and an empyt
all_min  = np.zeros(no_tests-1)

starting_point = 1

# now, let's run the tests!
for i in range(no_tests-1):
    if control_files[i] == "gsl":
        package = "gsl"
    else: # assume ral_nlls is being called
        package = "ral_nlls"
    
    try:
        subprocess.call(["cp", "control_files/"+control_files[i], \
                         cutestdir+"/src/"+package+"/"+package.upper()+".SPC"])
    except:
        print "Error: No control file " + control_files[i] + "found"
    
    os.chdir("cutest/sif/")

    subprocess.call(["cp","/dev/null",control_files[i]+"_iter.out"])
    if i == 0:
        # very first call, so create blank file...
        subprocess.call(["runcutest","-p",package,"--decode",problem,"-st",str(starting_point)])
    else: # no need to decode again....
        subprocess.call(["runcutest","-p",package])
    subprocess.call(["mv", control_files[i]+"_iter.out", \
                     "../../data/"+control_files[i]+"_iter.out"])

    os.chdir("../../")
    filename = "data/" + control_files[i] + "_iter.out"
    progress[i] = np.loadtxt(filename,dtype = datatype)
    all_min[i] = progress[i]['res'].min()

print all_min

minvalue = all_min.min()
mineps = minvalue - np.finfo(float).eps

#
# do the plotting!
#

plt.figure(1)
for i in range(no_tests-1):
    plt.semilogy(progress[i]['grad'], label=control_files[i])
plt.legend()

plt.title(problem+': gradients')
plt.xlabel('Iteration number')
plt.ylabel('$||J^Tr||_2$')
plt.savefig('../../../doc/img/'+problem+'_'+short_hash+'.png')

plt.figure(2)
for i in range(no_tests-1):
    plt.semilogy(progress[i]['res']-mineps, label=control_files[i])

plt.legend()
plt.title(problem+': residuals \n minimizer = '+str(minvalue) )
plt.xlabel('Iteration number')
plt.ylabel('$1/2||r_k||^2_2 - 1/2||r_*||^2_2$')
plt.savefig('../../../doc/img/'+problem+'_res_'+short_hash+'.png')

plt.show()

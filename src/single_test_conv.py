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
    
    os.chdir("../cutest/sif/")

    if i == 0:
        # very first call, so create blank file...
        subprocess.call(["cp","/dev/null",control_files[i]+"_iter.out"])
        subprocess.call(["runcutest","-p",package,"--decode",problem])
    else: # no need to decode again....
        subprocess.call(["runcutest","-p",package])
    subprocess.call(["mv", control_files[i]+"_iter.out", \
                     "../../src/data/"+control_files[i]+"_iter.out"])

    os.chdir("../../src/")
    filename = "data/" + control_files[i] + "_iter.out"
    progress[i] = np.loadtxt(filename,dtype = datatype)
    all_min[i] = progress[i]['res'].min()

print all_min

minvalue = all_min.min()
mineps = minvalue - np.finfo(float).eps

#
# do the plotting!
#

# plt.figure(1)
for i in range(no_tests-1):
    plt.figure(i)
    plt.semilogy(progress[i]['grad'], label='gradient')
    function_value = 0.5*((progress[i]['res'])**2)
    plt.semilogy(function_value, label='residual')
    plt.semilogy(progress[i]['grad']/function_value, label='ratio')
    plt.legend()
    plt.title(problem+': '+control_files[i])
    plt.xlabel('Iteration number')
plt.show()

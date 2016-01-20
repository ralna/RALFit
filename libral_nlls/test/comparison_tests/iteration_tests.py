#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os

# Let's get the files containing the problem and control parameters from the calling command..
no_tests = len(sys.argv)-1

print "*************************************"
print "**                                 **"
print "**        R A L _ N L L S          **"
print "**                                 **"
print "*************************************"

control_files = [None for i in range(no_tests)]
for i in range(no_tests):
    control_files[i] = sys.argv[i+1]

# read the cutest directory from the environment variables
try:
    cutestdir = os.environ.get("CUTEST")
except:
    print "ERROR: the CUTEST environment variable doesn't appear to be set"

# get the current git hash
short_hash = subprocess.check_output(['git','rev-parse','--short','HEAD']).strip()

# setup the datatype that we'll store the results in
info = np.dtype({'names' :   ['pname','n','m','status','iter',
                              'func','jac','hess',
                              'res','grad','ratio'],
                 'formats' : ['S10' ,int ,int,int,int,
                              int, int, int,
                              float,float,float]})

# get the list of problems...
#prob_list = "nist"
prob_list = "names_nist_first"
problems = np.loadtxt("cutest/sif/"+prob_list+".txt", dtype = str)

no_probs = len(problems)

for i in range(no_probs):
    print "**** "+ problems[i] +" ****"
    # now, let's run the tests!
    for j in range(no_tests):
        if control_files[j] == "gsl":
            package = "gsl"
        else: # assume ral_nlls is being called
            package = "ral_nlls"
            
        try:
            subprocess.call(["cp", "control_files/"+control_files[j], \
                             cutestdir+"/src/"+package+"/"+package.upper()+".SPC"])
        except:
            print "Error: No control file " + control_files[i] + "found"
           
        os.chdir("cutest/sif/")
               
        if i == 0:
            # very first call, so create blank file...
            subprocess.call(["cp","/dev/null",control_files[j]+".out"])

        if j == 0:
            # and then call sifdecoder as well as cutest
            subprocess.call(["runcutest","-p",package,"--decode",problems[i]])
        else: # no need to decode again....
            subprocess.call(["runcutest","-p",package])
        
        os.chdir("../../")

# now we have all the data, we just need to process it....

data = [None for i in range(no_tests)]
clear_best = np.zeros(no_tests, dtype = np.int)
best = np.zeros(no_tests,dtype = np.int)
too_many_its = np.zeros(no_tests, dtype = np.int)

for j in range(no_tests):
    subprocess.call(["mv", "cutest/sif/"+control_files[j]+".out", \
                     "data/"+control_files[j]+".out"])
    data[j] = np.loadtxt("data/"+control_files[j]+".out", dtype = info)
    if control_files[j] == "gsl":
        too_many_its[j] = -2
    else:
        too_many_its[j] = -1

all_iterates = [data[j]['iter'] for j in range(no_tests)]
all_status = [data[j]['status'] for j in range(no_tests)]

local_iterates = np.zeros(no_tests, dtype = np.int)

# finally, run through the data....
for i in range(0,no_probs):
    for j in range (0,no_tests):
        if (all_status[j][i] != 0) and (all_status[j][i] != too_many_its[j]):
            all_iterates[j][i] = 9999
        local_iterates[j] = all_iterates[j][i]
    minvalue = local_iterates.min()
    if (minvalue == 9999) or (minvalue == 1000): continue
    minima = np.where( local_iterates == minvalue )
    if minima[0].shape[0] == 1:
        clear_best[ minima[0][0] ] += 1
    for j in range(0,minima[0].shape[0]):
        best[ minima[0][j] ] += 1

print "Iteration numbers, git commit "+short_hash
print "%10s" % "problem",
for j in range(0,no_tests):
    print " ", 
    print "%16s" % control_files[j],
print " "
    
for i in range(0,no_probs):
    print "%10s" % problems[i],
    for j in range(0,no_tests):
        print ' ', 
        print "%16d" % all_iterates[j][i],
    print ' '

print "\n\n"
for j in range (0, no_tests):
    print control_files[j]+" is best ",best[j],\
        " times (and clear best ",clear_best[j]," times)"
    
    

#                    filename = "data/" + control_files[i] + "_iter.out"
#                    progress[i] = np.loadtxt(filename,dtype = datatype)
#                    all_min[i] = progress[i]['res'].min()
                    
#print all_min
#
#minvalue = all_min.min()
#mineps = minvalue - np.finfo(float).eps

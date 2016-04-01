#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import analyse_data

def analyse(no_tests,no_probs,control_files):
    # get the current git hash
    short_hash = subprocess.check_output(['git','rev-parse','--short','HEAD']).strip()

    # setup the datatype that we'll store the results in
    info = np.dtype({'names' :   ['pname','n','m','status','iter',
                                  'func','jac','hess',
                                  'res','grad','ratio'],
                     'formats' : ['S10' ,int ,int,int,int,
                                  int, int, int,
                                  float,float,float]})

    data = [None for i in range(no_tests)]
    clear_best = np.zeros(no_tests, dtype = np.int)
    best = np.zeros(no_tests,dtype = np.int)
    too_many_its = np.zeros(no_tests, dtype = np.int)
    local_iterates = np.zeros(no_tests, dtype = np.int)
    average_iterates = np.zeros(no_tests, dtype = np.int)
    average_funeval = np.zeros(no_tests, dtype = np.int)
    no_failures = np.zeros(no_tests, dtype = np.int)

    for j in range(no_tests):
        subprocess.call(["mv", "cutest/sif/"+control_files[j]+".out", \
                         "data/"+control_files[j]+".out"])
        data[j] = np.loadtxt("data/"+control_files[j]+".out", dtype = info)
        if control_files[j] == "gsl":
            too_many_its[j] = -2
        else:
            too_many_its[j] = -1

    all_iterates = [data[j]['iter'] for j in range(no_tests)]
    all_func = [data[j]['func'] for j in range(no_tests)]
    all_status = [data[j]['status'] for j in range(no_tests)]

    # finally, run through the data....
    for i in range(0,no_probs):
        for j in range (0,no_tests):
            if (all_status[j][i] != 0) and (all_status[j][i] != too_many_its[j]):
                all_iterates[j][i] = -9999 
            local_iterates[j] = all_iterates[j][i]
            if (all_iterates[j][i] < 0):
                no_failures[j] += 1
            else:
                average_iterates[j] += all_iterates[j][i]
                average_funeval[j] += all_func[j][i]
        minvalue = np.absolute(local_iterates).min()
        if (minvalue == 9999) or (minvalue == 1000): continue
        minima = np.where( local_iterates == minvalue )
        if minima[0].shape[0] == 1:
            clear_best[ minima[0][0] ] += 1
        for j in range(0,minima[0].shape[0]):
            best[ minima[0][j] ] += 1

    for j in range(0,no_tests):
        average_funeval[j] = average_funeval[j] / (no_probs - no_failures[j])
        average_iterates[j] = average_iterates[j] / (no_probs - no_failures[j])

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

    for j in range (0, no_tests):
        print control_files[j]+" took ",average_iterates[j],\
            " iterations and ", average_funeval[j], \
            " func. evals on average, and failed ", no_failures[j]," times)"
        return;

if __name__ == "__main__":
    analyse()

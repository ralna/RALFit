#!/usr/bin/python

import numpy as np
import subprocess

short_hash =  subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])

info = np.dtype({'names' :   ['pname','n','m','status','iter','res','grad'],
                 'formats' : ['S10' ,int ,int,int,     int,   float,float]})

data1 = np.loadtxt("data/results_Deltabase.out", dtype = info)
data2 = np.loadtxt("data/results_Delta-3.out", dtype = info)
data3 = np.loadtxt("data/results_Delta0.out", dtype = info)
data4 = np.loadtxt("data/results_Delta3.out", dtype = info)


#problems = np.loadtxt("sif_names.txt", dtype=str)
problems = np.array(data1['pname']);

all_iterates = np.array([ data1['iter'], 
                          data2['iter'], 
                          data3['iter'],
                          data4['iter'] ])
all_status = np.array([ data1['status'], 
                        data2['status'], 
                        data3['status'],
                        data4['status'] ])

clear_best = [0, 0, 0, 0]
best = [0, 0, 0, 0]
too_many_its = [-1, -1, -1, -1]

for i in range(0,data1.shape[0]):
    for j in range(0,4):
        # set the ones the failed to 9999 so they are not counted
        if (all_status[j,i] != 0) and (all_status[j,i] != too_many_its[j] ):
            all_iterates[j,i] = 9999

    minvalue = all_iterates[:,i].min()
    if (minvalue == 9999) or (minvalue == 1000): continue
    minima = np.where( all_iterates[:,i] == minvalue )
    if minima[0].shape[0] == 1:
        clear_best[ minima[0][0] ] += 1
    for j in range(0,minima[0].shape[0]):
        best[ minima[0][j] ] += 1


for i in range(0,data1.shape[0]):
    print problems[i],' ',all_iterates[0,i],' ', all_iterates[1,i],' ', \
          all_iterates[2,i],' ', all_iterates[3,i]

print "\n\nBase is best ",best[0]," times (and clear best ",clear_best[0]," times)"
print "1.0E-3 is best ",best[1]," times (and clear best ",clear_best[1]," times)"
print "1.0 is best ",best[2]," times (and clear best ",clear_best[2]," times)"
print "1.0E3 is best ",best[3]," times (and clear best ",clear_best[3]," times)"


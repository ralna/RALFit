#!/usr/bin/python

import numpy as np
import subprocess

short_hash =  subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])

info = np.dtype({'names' :   ['pname','n','m','status','iter','res','grad'],
                 'formats' : ['S10' ,int ,int,int,     int,   float,float]})

data1 = np.loadtxt("data/results_m1.out", dtype = info)
data2 = np.loadtxt("data/results_m2.out", dtype = info)
data7 = np.loadtxt("data/results_m7.out", dtype = info)
data8 = np.loadtxt("data/results_m8.out", dtype = info)

#problems = np.loadtxt("sif_names.txt", dtype=str)
problems = np.array(data1['pname']);

all_iterates = np.array([ data1['iter'], 
                          data2['iter'], 
                          data7['iter'],
                          data8['iter'] ])
all_status = np.array([ data1['status'], 
                        data2['status'], 
                        data7['status'],
                        data8['status'],])

best = [0, 0, 0, 0]

for i in range(0,data1.shape[0]):
    for j in range(0,4):
        # set the ones the failed to 9999 so they are not counted
        if (all_status[j,i] != 0) and (all_status[j,i] != -1):
            all_iterates[j,i] = 9999
    minvalue = all_iterates[:,i].min()
    if (minvalue == 9999) or (minvalue == 1000): continue
    minima = np.where( all_iterates[:,i] == minvalue )
    for j in range(0,minima[0].shape[0]):
        best[ minima[0][j] ] += 1


for i in range(0,data1.shape[0]):
    print problems[i],' ',all_iterates[0,i],' ', all_iterates[1,i],' ', \
          all_iterates[2,i],' ', all_iterates[3,i]  

print "\n\n Gauss-Newton is best ",best[0]," times"
print "Newton is best ",best[1]," times"
print "Hybrid is best ",best[2]," times"
print "Hybrid II is best ",best[3]," times"

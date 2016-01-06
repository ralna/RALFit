#!/usr/bin/python

import numpy as np

info = np.dtype({'names' : ('n','m','status','iter','res','grad'),
        'formats' : ('int','int','int','int','float','float')})

data1 = np.loadtxt("data/results_m1.out_data", dtype = info)
data2 = np.loadtxt("data/results_m2.out_data", dtype = info)
data7 = np.loadtxt("data/results_m7.out_data", dtype = info)

all_iterates = np.array([ data1['iter'] , data2['iter'], data7['iter']  ])
all_status = np.array([ data1['status'] , data2['status'], data7['status']  ])

best = [0, 0, 0]

for i in range(0,data1.shape[0]):
    for j in range(0,3):
        # set the ones the failed to 9999 so they are not counted
        if all_status[j,i] != 0:
            all_iterates[j,i] = 9999
    minvalue = all_iterates[:,i].min()
    if minvalue == 9999: continue
    minima = np.where( all_iterates[:,i] == minvalue )
    for j in range(0,minima[0].shape[0]):
        best[ minima[0][j] ] += 1

print "Gauss-Newton is best ",best[0]," times"
print "Newton is best ",best[1]," times"
print "Hybrid is best ",best[2]," times"

for i in range(0,data1.shape[0]):
    print all_iterates[0,i],' ', all_iterates[1,i],' ', all_iterates[2,i]

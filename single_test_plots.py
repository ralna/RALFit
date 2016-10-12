#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import argparse

def main():
    # Let's get the files containing the problem and control parameters from the calling command..
    parser = argparse.ArgumentParser()
    parser.add_argument('control_files', nargs='*', help="the names of lists of optional arguments to be passed to nlls_solve, in the format required by CUTEST, which are found in files in the directory ./control_files/")
    parser.add_argument("-p","--problem", help="which problem should be use?")
    parser.add_argument("-n","--no_plot", help="if present, we do not draw the plots", action="store_true")
    parser.add_argument("-s","--starting_point", help="which starting point to use? (default=1)",default=1)
    args = parser.parse_args()
    control_files = args.control_files
    no_tests = len(control_files) #!len(sys.argv)-1
    problem = args.problem

    print "Solving problem " + problem + " (starting point "+str(args.starting_point)+")"

    # read the cutest directory from the environment variables
    try:
        cutestdir = os.environ.get("CUTEST")
    except:
        print "ERROR: the CUTEST environment variable doesn't appear to be set"
    
    # copy the ral_nlls files to cutest
    subprocess.call(["cp","cutest/src/ral_nlls/ral_nlls_test.f90",cutestdir+"/src/ral_nlls/"])
    subprocess.call(["cp","cutest/src/ral_nlls/ral_nlls_main.f90",cutestdir+"/src/ral_nlls/"])
        
    # get the current git hash
    short_hash = subprocess.check_output(['git','rev-parse','--short','HEAD']).strip()

    # setup the datatype that we'll store the results in
    datatype = np.dtype({'names' : ('res','grad'),
                         'formats' : ('float','float')})

    # setup a python list in which to store the arrays...
    progress = [None for i in range(no_tests)]
    # and an empyt
    all_min  = np.zeros(no_tests)

    # now, let's run the tests!
    for i in range(no_tests):
        if control_files[i] == "gsl":
            package = "gsl"
        else: # assume ral_nlls is being called
            package = "ral_nlls"

        try:
            subprocess.call(["cp", "control_files/"+control_files[i], \
                             "cutest/sif/"+package.upper()+".SPC"])
        except:
            print "Error: No control file " + control_files[i] + "found"

        os.chdir("cutest/sif/")

        subprocess.call(["cp","/dev/null",control_files[i]+"_iter.out"])
        if i == 0:
            # very first call, so create blank file...
            subprocess.call(["runcutest","-p",package,"--decode",problem,"-st",str(args.starting_point)])
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
    if not args.no_plot:
        plot(no_tests,control_files,progress,short_hash,problem,mineps,minvalue)

def plot(no_tests,control_files,progress,short_hash,problem,mineps,minvalue):
    plt.figure(1)
    for i in range(no_tests):
        plt.semilogy(progress[i]['grad'], label=control_files[i])
    plt.legend()

    print os.getcwd()
    plt.title(problem+': gradients')
    plt.xlabel('Iteration number')
    plt.ylabel('$||J^Tr||_2$')
    plt.savefig('data/img/'+problem+'_'+short_hash+'.png')

    plt.figure(2)
    for i in range(no_tests):
        plt.semilogy(progress[i]['res']-mineps, label=control_files[i])

    plt.legend()
    plt.title(problem+': residuals \n minimizer = '+str(minvalue) )
    plt.xlabel('Iteration number')
    plt.ylabel('$1/2||r_k||^2_2 - 1/2||r_*||^2_2$')
    plt.savefig('data/img/'+problem+'_res_'+short_hash+'.png')

    plt.show()


if __name__ == "__main__":
    main()

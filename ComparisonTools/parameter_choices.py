#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os


def main(argv):
    
    # we want to find good choices for
    # initial TR radius
    # eta_successful
    # eta_success_but_reduce
    # eta_very_successful
    # eta_too_successful
    # radius_increase
    # radius_reduce
    # hybrid_tol
    # hybrid_switch_its
    # (scale?)
    
    print "++ ral_nlls :: parameter optimization ++"
    
    initial_tr_radius_line = 10
    summary_file_line = 30
    iteration_summary_line = 32
    
    method_values = np.array([3,4])
    method_line = 4
    initial_radius_values = np.array([90.0, 100.0]) 
    initial_radius_line = 10
    scale_trim_min_values = np.array(["T"])
    scale_trim_min_line = 24
    scale_trim_max_values = np.array(["T"])
    scale_trim_max_line = 25
    scale_require_increase_values = np.array(["F"])
    scale_require_increase_line = 26

    
    # start at DEF
    control_files = []
    description = []
    no_tests = 0
    testnumber = 0
    for method in method_values:
        for initial_radius in initial_radius_values:
            for scale_trim_min in scale_trim_min_values:
                for scale_trim_max in scale_trim_max_values:
                    for scale_require_increase in scale_require_increase_values:
                        testnumber += 1
                        linenumber = 0
                        filename = "TEST"+str(testnumber)
                        desc_string = filename+":"
                        no_tests += 1
                        control_files.append(filename)
                        newfile = open("control_files/"+filename,"w")
                        with open("control_files/DEF","r") as f:
                            for line in f:
                                linenumber += 1
                                if linenumber == summary_file_line:
                                    newfile.write(filename+".out    \n")
                                elif linenumber == iteration_summary_line:
                                    newfile.write(filename+"_iter.out   ")
                                elif linenumber == method_line:
                                    newfile.write(str(method)+"\n")
                                elif linenumber == initial_radius_line:
                                    newfile.write(str(initial_radius)+"\n")
                                    desc_string = desc_string+" initial_radius="+str(initial_radius)
                                elif linenumber == scale_trim_max_line:
                                    newfile.write(str(scale_trim_max)+"\n")
                                    desc_string = desc_string+" scale_trim_max="+str(scale_trim_max)
                                elif linenumber == scale_trim_min_line:
                                    newfile.write(str(scale_trim_min)+"\n")
                                    desc_string = desc_string+" scale_trim_min="+str(scale_trim_min)
                                elif linenumber == scale_require_increase_line:
                                    newfile.write(str(scale_require_increase)+"\n")
                                    desc_string = desc_string+\
                                                  "scale_require_increase="+str(scale_require_increase)
                                else:
                                    newfile.write(line)
                        newfile.close()
                        description.append(desc_string)

            
    # get the problems that we'll be running....
    prob_list = "all_not_diamond"
    problems = np.loadtxt("cutest/sif/"+prob_list+".txt", dtype = str)
    no_probs = len(problems)

    for probindex, problem in enumerate(problems, start=0):
        print "*** "+ problem +" ***"
        compute(no_tests,control_files,problem,probindex)
    # this prints results onto the file control_file.out
    # Let's process it...
    
    info = np.dtype({'names' :   ['pname','n','m','status','iter',
                                  'func','jac','hess',
                                  'res','grad','ratio'],
                     'formats' : ['S10' ,int ,int,int,int,
                                  int, int, int,
                                  float,float,float]})
    hashinfo = np.dtype({'names'   : ['hash','no_probs'], 
                         'formats' : ['S7',int]})

    data = [None for i in range(no_tests)]
    metadata = [None for i in range(no_tests)]
    clear_best = np.zeros(no_tests, dtype = np.int)
    best = np.zeros(no_tests,dtype = np.int)
    too_many_its = np.zeros(no_tests, dtype = np.int)
    local_iterates = np.zeros(no_tests, dtype = np.int)
    average_iterates = np.zeros(no_tests, dtype = np.int)
    average_funeval = np.zeros(no_tests, dtype = np.int)
    no_failures = np.zeros(no_tests, dtype = np.int)
    no_dnc = np.zeros(no_tests, dtype = np.int)


    for index, control_file in enumerate(control_files, start=0):
        subprocess.call(["mv", "cutest/sif/"+control_file+".out", \
                             "data/"+control_file+".out"])
        # get the hash of the git version
        short_hash = subprocess.check_output(['git','rev-parse','--short','HEAD']).strip()
        f = open('data/'+control_file+".hash",'w')
        f.write(short_hash+"\t"+str(no_probs))
        f.close()
        data[index] = np.loadtxt("data/"+control_file+".out", dtype = info)
        metadata[index] = np.loadtxt("data/"+control_file+".hash", dtype = hashinfo)
        if control_file == "gsl":
            too_many_its[index] = -2
        else:
            too_many_its[index] = -1

    all_iterates = [data[j]['iter'] for j in range(no_tests)]
    all_func = [data[j]['func'] for j in range(no_tests)]
    all_status = [data[j]['status'] for j in range(no_tests)]

    # finally, run through the data....
    for j in range (0,no_tests):
        if j == 0:
            short_hash = str(metadata[j]['hash'])
            hash_error = False
        elif str(metadata[j]['hash']) != short_hash:
            hash_error = True
    
    for i in range(0,no_probs):     
        for j in range (0,no_tests):
            if (all_status[j][i] != 0) and (all_status[j][i] != too_many_its[j]):
                all_iterates[j][i] = 9999 
                no_failures[j] += 1
            elif (all_status[j][i] == too_many_its[j]):
                all_iterates[j][i] = - all_iterates[j][i]   
                all_func[j][i] = - all_func[j][i]
                no_dnc[j] += 1 
            local_iterates[j] = all_iterates[j][i]
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
        average_funeval[j] = average_funeval[j] / no_probs
        average_iterates[j] = average_iterates[j] / no_probs


    for i in range(0,no_probs):
        print "%10s" % problems[i],
        for j in range(0,no_tests):
            print ' ', 
            print "%16d" % all_iterates[j][i],
        print ' '

    print "\n\n"
    for j in range (0, no_tests):
        print description[j]+" is best ",best[j],\
            " times (and clear best ",clear_best[j]," times)"

    for j in range (0, no_tests):
        print description[j]+" took ",average_iterates[j],\
            " iterations and ", average_funeval[j], \
            " func. evals on average. Method failed ", no_failures[j]," times", \
            " and did not converge ", no_dnc[j]," times"
        
        
# ++ define the compute function ++ #

def compute(no_tests,control_files,problem,probindex):
    # read the cutest directory from the environment variables
    try:
        cutestdir = os.environ.get("CUTEST")
    except:
        raise Error("the CUTEST environment variable doesn't appear to be set")
    
    starting_point = 1
    
    for index, control_file in enumerate(control_files, start=0):
        print "*** control file = "+control_file+" ***"
        if control_file == "gsl":
            package = "gsl"
        else: # assume ral_nlls is being called
            package = "ral_nlls"
            
        try:
            subprocess.call(["cp", "control_files/"+control_file, \
                             "cutest/sif/"+package.upper()+".SPC"])
        except:
            print "The control file (control_files/"+control_file+") is not found"
            raise 

           
        os.chdir("cutest/sif/")
               
        if probindex == 0:
            # very first call, so create blank file...
            subprocess.call(["cp","/dev/null",control_file+".out"])

        if index == 0:
            # and then call sifdecoder as well as cutest
            subprocess.call(["runcutest","-p",package,"--decode",problem, \
                             "-st",str(starting_point)])
        else: # no need to decode again....
            subprocess.call(["runcutest","-p",package])
        
        os.chdir("../../")

if __name__ == "__main__":
    main(sys.argv[1:])

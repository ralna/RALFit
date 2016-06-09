#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os

def main(argv):
    # Let's get the files containing the problem and control parameters from the calling command..
    no_tests = len(argv)
    if argv[-1] == 'no':
        no_tests = no_tests - 1
        compute_results = False
    else:
        compute_results = True

    print "*************************************"
    print "**                                 **"
    print "**        R A L _ N L L S          **"
    print "**                                 **"
    print "*************************************"

    # put the names of the control files to be run into an array
    control_files = [None for i in range(no_tests)]
    for i in range(no_tests):
        control_files[i] = argv[i]

    # get the list of problems...
#    prob_list = "nist"
#    prob_list = "nist_average"
    prob_list = "names_nist_first"
    #prob_list = "names_minus_boxbod"
    problems = np.loadtxt("cutest/sif/"+prob_list+".txt", dtype = str)

    no_probs = len(problems)

    print "Testing ",no_probs," problems with ",no_tests," minimizers"

    for i in range(no_probs):
        if compute_results:
            # let's run the tests!
            print "**** "+ problems[i] +" ****"
            compute(no_tests,control_files,problems,i)
    
    if compute_results:
        for j in range(no_tests):
            subprocess.call(["mv", "cutest/sif/"+control_files[j]+".out", \
                             "data/"+control_files[j]+".out"])
            if control_files[j] == "gsl":
                # add a newline
                gslresults = open("data/gsl.out","a")
                gslresults.write("\n")
                gslresults.close()
            # get the hash of the git version
            short_hash = subprocess.check_output(['git','rev-parse','--short','HEAD']).strip()
            f = open('data/'+control_files[j]+".hash",'w')
            f.write(short_hash+"\t"+str(no_probs))
            f.close()

    # # now we have all the data, we just need to process it....
     
    # setup the datatype that we'll store the results in
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

    for j in range(no_tests):
        data[j] = np.loadtxt("data/"+control_files[j]+".out", dtype = info)
        metadata[j] = np.loadtxt("data/"+control_files[j]+".hash", dtype = hashinfo)
        if control_files[j] == "gsl":
            too_many_its[j] = -2
        else:
            too_many_its[j] = -1

    all_iterates = [data[j]['iter'] for j in range(no_tests)]
    all_func = [data[j]['func'] for j in range(no_tests)]
    all_status = [data[j]['status'] for j in range(no_tests)]
    
    normalized_mins = [data[j]['res'] for j in range(no_tests)]
    tiny = 1e-8
    normalized_iterates = np.copy(all_iterates)#[data[j]['iter'] for j in range(no_tests)]
    normalized_func = np.copy(all_func)
    failure = np.zeros((no_probs, no_tests))

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
                all_iterates[j][i] = -9999 
                failure[i][j] = 1 
                normalized_iterates[j][i] = 9999
                normalized_func[j][i] = 9999
            local_iterates[j] = all_iterates[j][i]
            if (all_iterates[j][i] < 0):
                no_failures[j] += 1
                if (failure[i][j] != 1):
                    failure[i][j] = 2
                    normalized_iterates[j][i] = -normalized_iterates[j][i]
                    normalized_func[j][i] = -normalized_func[j][i]
            else:
                average_iterates[j] += all_iterates[j][i]
                average_funeval[j] += all_func[j][i]
            if normalized_mins[j][i] < tiny:
                # truncate anything smaller than tiny
                normalized_mins[j][i] = tiny 
#            normalized_mins[j][i] = normalized_mins[j][i]/smallest_resid[i]
#            normalized_iterates[j][i] = normalized_iterates[j][i]  + 1 - smallest_iterates[i]
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

    smallest_resid = np.amin(normalized_mins, axis = 0)
    smallest_iterates = np.amin(np.absolute(normalized_iterates), axis = 0)        
    smallest_func = np.amin(np.absolute(normalized_func), axis = 0)        
    normalized_mins = np.transpose(normalized_mins)
    normalized_iterates = np.transpose(normalized_iterates)
    normalized_func = np.transpose(normalized_func)
    mins_boundaries = np.array([1.1, 1.33, 1.75, 3.0])
    iter_boundaries = np.array([2, 5, 10, 30])
    func_boundaries = np.array([2, 5, 10, 30])
    additive = 0
    print_to_html(no_probs, no_tests, problems, normalized_mins, smallest_resid, 
                  mins_boundaries, 'normalized_mins', control_files, failure, additive)
    additive = 1
    print_to_html(no_probs, no_tests, problems, normalized_iterates, smallest_iterates, 
                  iter_boundaries, 'normalized_iters', control_files, failure, additive)
    print_to_html(no_probs, no_tests, problems, normalized_func, smallest_func, 
                  func_boundaries, 'normalized_func', control_files, failure, additive)
    
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

    if hash_error == True:
        print "\n\n"
        print "************************************************"
        print "*                 W A R N I N G               **"
        print "* results computed with different git commits  *"
        print "************************************************"
        print "\n"

    plot_prof(control_files,no_tests,prob_list)

def print_to_html(no_probs, no_tests, problems, data, smallest, boundaries, 
                  filename, control_files, failure, additive):

    # first, let's set the background colours...
    good = '#00ff00'
    averagegood = '#7fff00'
    average = '#ffff00'
    badaverage = '#ff7f00'
    bad = '#ff0000'

    print filename
    output = open('data/'+filename+'.html','w')
    output.write('<!DOCTYPE html>\n')
    output.write('<html>\n')
    output.write('<body>\n \n')

    output.write('Data is shown below.  The best value is shown in ')
    output.write('<span style=background-color:'+good+'>green</span> ')
    output.write('and the others are colour-coded depending on the size of ')
    if additive:
        output.write('|value - best|')
    else:
        output.write('|value / best|')
    output.write(', as shown in the key:\n')
    output.write('<table>\n')
    output.write('  <tr>\n')
    output.write('    <td bgcolor = '+good+'> x &lt; '+str(boundaries[0])+'</td>\n')
    output.write('    <td bgcolor = '+averagegood+'>') 
    output.write(str(boundaries[0])+' &le; x &lt; '+str(boundaries[1])+' </td>\n')
    output.write('    <td bgcolor = '+average+'>') 
    output.write(str(boundaries[1])+' &le; x &lt; '+str(boundaries[2])+' </td>\n')
    output.write('    <td bgcolor = '+badaverage+'>')
    output.write(str(boundaries[2])+' &le; x &lt; '+str(boundaries[3])+' </td>\n')
    output.write('    <td bgcolor = '+bad+'>  x &ge; '+str(boundaries[3])+' </td>\n')
    output.write('  </tr>\n')
    output.write('</table>\n')

    output.write('<table>\n')
    output.write('  <tr>\n')
    output.write('    <td>&dagger; denotes problems where the method failed </td>\n')
    output.write('  </tr>\n')
    output.write('  <tr>\n')
    output.write('    <td>&Dagger; denotes problems where the max number of iterations was reached </td> \n')
    output.write('  </tr>\n')
    output.write('</table>\n')
    
    output.write('<table>\n')
    output.write('  <tr>\n')
    output.write('    <td></td>\n')
    for j in range(0,no_tests):
        output.write('    <td>'+control_files[j]+'</td>\n')
    output.write('  </tr>\n')
        

    for i in range(0,no_probs):
        output.write('  <tr>\n')
        output.write('    <td>'+problems[i]+'</td>')
        for j in range(0,no_tests):
            if additive:
                current_value = data[i][j] - smallest[i] + 1
            else:
                current_value = data[i][j] / smallest[i]
            if failure[i][j] == 1:
                label = '&dagger;'
            elif failure[i][j] == 2:
                label = '&Dagger;'
            else:
                label = ''
            if current_value < boundaries[0]:
                colour = good
            elif current_value < boundaries[1]:
                colour = averagegood
            elif current_value < boundaries[2]:
                colour = average
            elif current_value < boundaries[3]:
                colour = badaverage
            else:
                colour = bad
            output.write('    <td bgcolor='+colour+'>'+str(format(data[i][j]))+label+'</td>')
        output.write('\n')
        output.write('  </tr>\n')
    output.write('</table>\n\n')
    output.write('</body>\n')
    output.write('</html>\n')
    output.close()

def format(number):
    if isinstance(number, int):
        return number
    elif isinstance(number, float):
        return "%.4g" % number



def compute(no_tests,control_files,problems,i):
    # read the cutest directory from the environment variables
    try:
        cutestdir = os.environ.get("CUTEST")
    except:
        raise Error("the CUTEST environment variable doesn't appear to be set")
    
    starting_point = 1

    for j in range(no_tests):
        if control_files[j] == "gsl":
            package = "gsl"
        else: # assume ral_nlls is being called
            package = "ral_nlls"
            
        try:
            subprocess.call(["cp", "control_files/"+control_files[j], \
                             "cutest/sif/"+package.upper()+".SPC"])
        except:
            raise Error("No control file " + control_files[j] + " found")
           
        os.chdir("cutest/sif/")
               
        if i == 0:
            # very first call, so create blank file...
            subprocess.call(["cp","/dev/null",control_files[j]+".out"])

        if j == 0:
            # and then call sifdecoder as well as cutest
            subprocess.call(["runcutest","-p",package,"--decode",problems[i], \
                             "-st",str(starting_point)])
        else: # no need to decode again....
            subprocess.call(["runcutest","-p",package])
        
        os.chdir("../../")

def plot_prof(control_files,no_tests,prob_list):
    # try:
    #     import pymatlab
    #     py_prof = True
    # except:
    #     print "If matlab is installed, install pymatlab\n"
    #     print " sudo pip install pymatlab\n"
    #     print "to allow performance profiles to be plotted natively\n"
    #     py_prof = False
    
    # performance profiles for iterations
    data_files = ""
    for j in range(no_tests):
        data_files += control_files[j]+".out"
        if j != no_tests-1:
            data_files += " "
    if prob_list=="names_nist_first" or prob_list=="sif_names":
        testset = "(All tests)"
    elif prob_list=="nist":
        testset = "NIST tests"
    else:
        testset = "CUTEst tests"

    os.chdir("data")
    try:
        subprocess.call(["pprof","5","iterations",data_files,testset])
        subprocess.call(["pprof","6","fevals",data_files,testset])
    except:
        print "Performance profiles not available: ensure pprof is in the path"

    os.chdir("..")

if __name__ == "__main__":
    main(sys.argv[1:])

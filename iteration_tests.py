#!/usr/bin/python
from __future__ import print_function, unicode_literals

import sys
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import argparse

def main():
    # Let's get the files containing the problem and control parameters from the calling command..
    parser = argparse.ArgumentParser()
    parser.add_argument('control_files', 
                        nargs='*', 
                        help="the names of lists of optional arguments to be passed to nlls_solve, in the format required by CUTEST, which are found in files in the directory ./control_files/")
    parser.add_argument("-r",
                        "--reuse_data", 
                        help="if present, we regenerate tables of iterations from previously computed data", 
                        action="store_true")
    parser.add_argument("-t",
                        "--test_times",
                        help="if present, runs tests ten times to produce accurate timing pprofs",
                        type=int,
                        default=0)
    parser.add_argument("-p",
                        "--problem_list", 
                        help="which list of problems to use? (default=names_nist_first)",
                        default="names_nist_first")
    parser.add_argument("-s",
                        "--starting_point", 
                        help="which starting point to use? (default=1)", 
                        default=1)
    parser.add_argument("-np",
                        "-no_performance_profile",
                        action="store_true",
                        help="if present, do not display the performance profile")
    args = parser.parse_args()
    no_tests = len(args.control_files)
    compute_results = not args.reuse_data
    
    prob_list = args.problem_list
    problems = np.genfromtxt("cutest/sif/"+prob_list+".txt", dtype = 'U10')
    no_probs = len(problems)

    print("*************************************")
    print("**                                 **")
    print("**        R A L _ N L L S          **")
    print("**                                 **")
    print("*************************************")

    print("Testing {} problems with {} minimizers".format(no_probs,no_tests))

    if compute_results:
        # run the tests!
        run_cutest_and_copy_results_locally(args,problems)

    # # now we have all the data, we just need to process it....
     
    # setup the datatype that we'll store the results in
    info = np.dtype({'names' :   ['pname','n','m','status','iter',
                                  'func','jac','hess','inner',
                                  'res','grad','ratio','solve_time'],
                     'formats' : ['U10' ,int ,int,int,int,
                                  int, int, int, int,
                                  float,float,float,float]})
    info_noinner = np.dtype({'names' :   ['pname','n','m','status','iter',
                                          'func','jac','hess',
                                          'res','grad','ratio','solve_time'],
                             'formats' : ['U10' ,int ,int,int,int,
                                          int, int, int, 
                                          float,float,float,float]})
    hashinfo = np.dtype({'names'   : ['hash','no_probs'], 
                         'formats' : ['U10',int]})

    data = [None for i in range(no_tests)]
    metadata = [None for i in range(no_tests)]
    clear_best = np.zeros(no_tests, dtype = np.int)
    best = np.zeros(no_tests,dtype = np.int)
    too_many_its = np.zeros(no_tests, dtype = np.int)
    local_iterates = np.zeros(no_tests, dtype = np.int)
    local_inner_it = np.zeros(no_tests,dtype = np.int)
    average_iterates = np.zeros(no_tests, dtype = np.int)
    average_funeval = np.zeros(no_tests, dtype = np.int)
    average_inner = np.zeros(no_tests,dtype = np.int)
    no_failures = np.zeros(no_tests, dtype = np.int)

    InnerResults = 1 # if gsl (or another method with no inners) is present, 
                     # then do not collect inner data (maybe fix)
    for j in range(no_tests):
        try:
            data[j] = np.genfromtxt("data/"+args.control_files[j]+".out", dtype = info)
        except ValueError:
            # these are results that don't include inner iterations
            # (i.e. from gsl)
            # only in this case, don't look for inner iterations
            data[j] = np.genfromtxt("data/"+args.control_files[j]+".out", dtype = info_noinner)
            InnerResults = 0
            # we want to put this back into a file that *does* have inner iterations,
            # so that the performance profiles work below...
            add_inner_information(args.control_files[j],data[j])
        metadata[j] = np.genfromtxt("data/"+args.control_files[j]+".hash", dtype = hashinfo)
        if "gsl" in args.control_files[j].lower(): # == "gsl":
            too_many_its[j] = -2
        else:
            too_many_its[j] = -1
   
        
    no_probs = len(data[0]['res'])

    if args.test_times > 0:
        # repeat the experiment another 'args.test_times' times, and get the average time over
        # that many runs
        local_data = [None for i in range(no_tests)]
        sum_times = np.zeros([no_tests, no_probs])
        for j in range(no_tests):
            sum_times[j][:] = data[j]['solve_time'][:]

        print("{}".format(sum_times))
        no_runs = 0
        while no_runs < args.test_times:
            run_cutest_and_copy_results_locally(args,problems)
            no_runs += 1
            print("run {}/{}".format(no_runs,args.test_times))

            for j in range(no_tests):
                try:
                    data[j] = np.loadtxt("data/"+args.control_files[j]+".out", dtype = info)
                except ValueError:
                    # these are results that don't include inner iterations
                    # (i.e. from gsl)
                    # only in this case, don't look for inner iterations
                    data[j] = np.loadtxt("data/"+args.control_files[j]+".out", dtype = info_noinner)
                    InnerResults = 0
                    # we want to put this back into a file that *does* have inner iterations,
                    # so that the performance profiles work below...
                    add_inner_information(args.control_files[j],data[j])
                sum_times[j][:] += data[j]['solve_time'][:]
                print("sum_times = {}".format(sum_times))
                print("solve_time = {}".format(data[j]['solve_time'][:]))
#        print "final solve time = "
        for j in range(no_tests):
            data[j]['solve_time'][:] = sum_times[j][:]/no_runs
            print("final solve time = {}".format(data[j]['solve_time'][:]))
        
            
            

            
    all_iterates = [data[j]['iter'] for j in range(no_tests)]
    all_func = [data[j]['func'] for j in range(no_tests)]
    all_status = [data[j]['status'] for j in range(no_tests)]
    all_solve_time = [data[j]['solve_time'] for j in range(no_tests)]
    if InnerResults:
        all_inner = [data[j]['inner'] for j in range(no_tests)]
    else:
        # since there's no inner iterations, set the number of inner iterations 
        # to equal the number inner iterations
        all_inner = all_iterates
    
    
    normalized_mins = [data[j]['res'] for j in range(no_tests)]
    tiny = 1e-8
    normalized_iterates = np.copy(all_iterates)#[data[j]['iter'] for j in range(no_tests)]
    normalized_func = np.copy(all_func)
    normalized_inner = np.copy(all_inner)
    normalized_solve_time = np.copy(all_solve_time)

    print("no_tests = {}".format(no_tests))
    print("no_probs = {}".format(no_probs))

    
    print("********************")
    for j in range(no_tests):
        for i in range(no_probs):
            print("{}".format(data[j]['pname'][i]))
    print("********************")
    
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
                normalized_inner[j][i] = 9999
            local_iterates[j] = all_iterates[j][i]
            if (all_iterates[j][i] < 0):
                no_failures[j] += 1
                if (failure[i][j] != 1):
                    failure[i][j] = 2
                    normalized_iterates[j][i] = -normalized_iterates[j][i]
                    normalized_func[j][i] = -normalized_func[j][i]
                    normalized_inner[j][i] = -normalized_inner[j][i]
            else:
                average_iterates[j] += all_iterates[j][i]
                average_funeval[j] += all_func[j][i]
                average_inner[j] += all_inner[j][i]
            if normalized_mins[j][i] < tiny:
                # truncate anything smaller than tiny
                normalized_mins[j][i] = tiny 
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
        average_inner[j] = average_inner[j] / (no_probs - no_failures[j])

    smallest_resid = np.amin(normalized_mins, axis = 0)
    smallest_iterates = np.amin(np.absolute(normalized_iterates), axis = 0)
    smallest_func = np.amin(np.absolute(normalized_func), axis = 0)
    smallest_inner = np.amin(np.absolute(normalized_inner), axis = 0)
    smallest_solve_time = np.amin(np.absolute(normalized_solve_time), axis = 0)
    normalized_mins = np.transpose(normalized_mins)
    normalized_iterates = np.transpose(normalized_iterates)
    normalized_func = np.transpose(normalized_func)
    normalized_inner = np.transpose(normalized_inner)
    normalized_solve_time = np.transpose(normalized_solve_time)
    mins_boundaries = np.array([1.1, 1.33, 1.75, 3.0])
    iter_boundaries = np.array([2, 5, 10, 30])
    func_boundaries = np.array([2, 5, 10, 30])
    inner_boundaries = np.array([2, 5, 10, 30])
    solve_time_boundaries = np.array([1.1, 1.33, 1.75, 3.0])
    additive = 0
    print_to_html(no_probs, no_tests, data[0]['pname'], normalized_mins, smallest_resid, 
                  mins_boundaries, 'normalized_mins', args.control_files, failure, additive,
                  short_hash)
    print_to_html(no_probs, no_tests, data[0]['pname'], normalized_solve_time, smallest_solve_time, 
                  solve_time_boundaries, 'normalized_solve_time', args.control_files, 
                  failure, additive, short_hash)
    additive = 1
    print_to_html(no_probs, no_tests, data[0]['pname'], normalized_iterates, smallest_iterates,
                  iter_boundaries, 'normalized_iters', args.control_files, failure, additive,
                  short_hash)
    print_to_html(no_probs, no_tests, data[0]['pname'], normalized_func, smallest_func, 
                  func_boundaries, 'normalized_func', args.control_files, failure, additive,
                  short_hash)
    print_to_html(no_probs, no_tests, data[0]['pname'], normalized_inner, smallest_inner, 
                  inner_boundaries, 'normalized_inner', args.control_files, 
                  failure, additive, short_hash)
    
    print("Iteration numbers, git commit {}".format(short_hash))
    print("problem   ")
    for j in range(0,no_tests):
        print(" {0: <16}".format(args.control_files[j]))
    print(" ")

    for i in range(0,no_probs):
        print("{0: <10}".format(data[0]['pname'][i]))
        for j in range(0,no_tests):
            print(' {0: <16} '.format(all_iterates[j][i]))

    print("\n\n")
    for j in range (0, no_tests):
        print("{} is best {} times (and clear best {} times)".format(
            args.control_files[j],best[j],clear_best[j]))

    for j in range (0, no_tests):
        print("{} took {} iterations and {} func. evals on average, and failed {} times)"
            .format(args.control_files[j],average_iterates[j],average_funeval[j],no_failures[j]))
        if average_inner[j] > average_iterates[j]:
            print("{} took inner iterations on average".format(
            args.control_files[j],average_inner[j]))

    for j in range (0,no_tests):
        print("{} took {}s to solve all problems".format(args.control_files[j],np.sum(all_solve_time[:][j])))

    if hash_error == True:
        print("\n\n")
        print("************************************************")
        print("*                 W A R N I N G               **")
        print("* results computed with different git commits  *")
        print("************************************************")
        print("\n")

    plot_prof(args.control_files,no_tests,prob_list,args.np)

def print_to_html(no_probs, no_tests, problems, data, smallest, boundaries, 
                  filename, control_files, failure, additive, short_hash):

    # first, let's set the background colours...
    good = '#00ff00'
    averagegood = '#7fff00'
    average = '#ffff00'
    badaverage = '#ff7f00'
    bad = '#ff0000'

    print(filename)
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
        output.write('    <td>{}</td>'.format(problems[i]))
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
            output_string = '    <td bgcolor={}>{}{}</td>'.format(colour,format(data[i][j]),label)
            output.write(output_string)
        output.write('\n')
        output.write('  </tr>\n')
    output.write('</table>\n')
    output.write('Computed using git commit <b>#'+short_hash+'</b>\n\n')
    output.write('</body>\n')
    output.write('</html>\n')
    output.close()
    print("Check output on file://"+os.getcwd()+"/data/"+filename+".html")

def format(number):
    if isinstance(number, int):
        return number
    elif isinstance(number, float):
        return "%.4g" % number


def add_inner_information(method_name,data):
    # take the file method_name.out (which does not have inner information),
    # and add it in.
    # This is needed so that the output file is the same as that for RALFit, and
    # therefore pypprof can compare them
    copy_file = open("data/"+method_name+".out","w")
    space = "   "
    for i in range(len(data['pname'])):
        line_to_print = (data['pname'][i]+
                         space+
                         str(data['n'][i])+
                         space+
                         str(data['m'][i])+
                         space+
                         str(data['status'][i])+
                         space+
                         str(data['iter'][i])+
                         space+
                         str(data['func'][i])+
                         space+
                         str(data['jac'][i])+
                         space+
                         "0"+
                         space+
                         str(data['hess'][i])+
                         space+
                         str(data['res'][i])+
                         space+
                         str(data['grad'][i])+
                         space+
                         str(data['ratio'][i])+
                         space+
                         str(data['solve_time'][i])+
                         "\n"
        )
        copy_file.write(line_to_print)
    copy_file.close()
    
def run_cutest_and_copy_results_locally(args,problems):
#    no_tests = len(args.control_files)
    no_tests = 0 
    no_probs = len(problems)
    for i in range(no_probs):
        # let's run the tests!
        print("**** "+ problems[i] +" ****")
        compute(no_tests,args.control_files,problems,i,args.starting_point)
    
    
    for test in args.control_files:#range(no_tests):
        # copy the data file locally
        subprocess.call(["mv", "cutest/sif/"+test+".out", \
                         "data/"+test+".out"])
        if "gsl" in test.lower(): # == "gsl":
            # c doesn't end with a newline, so add one
            gslresults = open("data/gsl.out","a")
            gslresults.write("\n")
            gslresults.close()
        # get the hash of the git version
        short_hash = subprocess.check_output(['git','rev-parse','--short','HEAD']).strip()
        f = open('data/'+test+".hash",'w')
        output_string = "{} \t{}".format(short_hash,no_probs)
        f.write(output_string)
        f.close()

    

def compute(no_tests,control_files,problems,i,starting_point):
    # read the cutest directory from the environment variables
    try:
        cutestdir = os.environ.get("CUTEST")
    except:
        raise Error("the CUTEST environment variable doesn't appear to be set")

    # copy the ral_nlls files to cutest
    subprocess.call(["cp","cutest/src/ral_nlls/ral_nlls_test.f90",cutestdir+"/src/ral_nlls/"])
    subprocess.call(["cp","cutest/src/ral_nlls/ral_nlls_main.f90",cutestdir+"/src/ral_nlls/"])

    j = 0
    
    for test in control_files:
        j += 1
        if "gsl" in test.lower(): # == "gsl":
            package = "gsl"
        else: # assume ral_nlls is being called
            package = "ral_nlls"
            
        try:
            subprocess.call(["cp", "control_files/"+test, \
                             "cutest/sif/"+package.upper()+".SPC"])
        except:
            raise Error("No control file " + test + " found")
           
        os.chdir("cutest/sif/")
               
        if i == 0:
            # very first call, so create blank file...
            subprocess.call(["cp","/dev/null",test+".out"])

        if (int(starting_point) == 0) and (problems[i] in open('nist.txt').read()):
            # we have a nist problem
            range_sp = ['1','2']
        else:
            range_sp = [str(starting_point)]

        for sp in range_sp:
            if (j == 1) or (int(starting_point)==0): # and then call sifdecoder as well as cutest
                subprocess.call(["runcutest","-p",package,"--decode",problems[i], \
                                 "-st",sp])
                no_tests += 1
            else: # no need to decode again....
                subprocess.call(["runcutest","-p",package,"-st",sp])
                no_tests += 1

        os.chdir("../../")

def plot_prof(control_files,no_tests,prob_list,np):
    # performance profiles for iterations
    Strings = ["python pypprof -ne -c 5 -s iterations ",
               "python pypprof -ne -c 6 -s fevals ",
               "python pypprof -ne -nf --log -c 13 -s time "]
    data_files = ""
    for j in range(no_tests):
        data_files += control_files[j]+".out"
        if j != no_tests-1:
            data_files += " "
    Strings[:] = [string + data_files for string in Strings]

    if prob_list=="names_nist_first" or prob_list=="sif_names":
        testset = "'test problems'"
    elif prob_list=="nist":
        testset = "'NIST tests'"
    else:
        testset = "'CUTEst tests'"

    Strings[:] = [string + " -t " + testset for string in Strings]

    if np:
        Strings[:] = [string + " -np" for string in Strings]

    os.chdir("data")
    try:
        for string in Strings:
            os.system(string)
    except:
        print("Performance profiles not available: ensure pprof is in the path")
    os.chdir("..")

if __name__ == "__main__":
    main()
#    main(sys.argv[1:])

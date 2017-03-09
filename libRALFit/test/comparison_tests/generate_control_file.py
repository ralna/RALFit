#!/usr/bin/python

import argparse

def main():
    ## Let's get the arguments from the user
    parser = argparse.ArgumentParser()
    parser.add_argument("control_name",
                        help="The name of the control file")
    parser.add_argument("--base_file","-b",
                        help="The base file to amend",default="DEF")
    opt_arguments = ['error',
                     'out',
                     'print_level',
                     'maxit',
                     'model',
                     'type_of_method',
                     'nlls_method',
                     'lls_solver',
                     'stop_g_absolute',
                     'stop_g_relative',
                     'stop_f_absolute',
                     'stop_f_relative',
                     'stop_s',
                     'relative_tr_radius',
                     'initial_radius_scale',
                     'initial_radius',
                     'base_regularization',
                     'regularization',
                     'regularization_term',
                     'regularization_power',
                     'maximum_radius',
                     'eta_successful',
                     'eta_success_but_reduce',
                     'eta_very_successful',
                     'eta_too_successful',
                     'radius_increase',
                     'radius_reduce',
                     'radius_reduce_max',
                     'tr_update_strategy',
                     'hybrid_switch',
                     'exact_second_derivatives',
                     'subproblem_eig_fact',
                     'use_ews_subproblem',
                     'scale',
                     'scale_max',
                     'scale_min',
                     'scale_trim_min',
                     'scale_trim_max',
                     'scale_require_increase',
                     'calculate_svd_J',
                     'setup_workspaces',
                     'remove_workspaces',
                     'more_sorensen_maxits',
                     'more_sorensen_shift',
                     'more_sorensen_tiny',
                     'more_sorensen_tol',
                     'hybrid_tol',
                     'hybrid_switch_its',
                     'reg_order',
                     'inner_method',
                     'output_progress_vectors',
                     'update_lower_order',
                     'summary_unit',
                     'iteration_summary']
    
    for arg in opt_arguments:
        parser.add_argument('--'+arg)
        
    options = parser.parse_args()
    options_dict = vars(options)

    print "++++ Generating a cutest control file for ral_nlls ++++"

    # Check we're not overwriting the defaults (which would be bad)
    if options.control_name == "DEF":
        print "** Don't overwrite DEF using this script **"
        exit()

    # read base_file in as a string
    options_string = open("control_files/"+options.base_file).read()
    
    string_length = 22           # this is the space we have before the control name

    # setup a string that will overwrite the default option
    blank = " "*string_length    
    filename_insert = options.control_name+blank[len(options.control_name):]

    new_values = 0

    for key in options_dict:
        if key == "base_file":
            continue # we don't need to write this option to file
        if options_dict[key] != None:
            # if the defaults have been changed, then 
            # edit the string accordingly
            control_location = options_string.find(key)
            keystr = str(options_dict[key])
            key_insert = keystr+blank[len(keystr):]
            old_value = options_string[control_location - string_length:control_location]
            print key+": replacing "+old_value.strip()+" with "+keystr
            if old_value != key_insert:
                new_values += 1
            options_string = options_string[:control_location - string_length]+ \
                             key_insert + \
                             options_string[control_location:]

    # check if the defaults have actually been changed
    if new_values < 2:
        print "You haven't changed the defaults at all"
        print "exiting without saving anything..."
        exit()

    # and, finally, write to a file
    file = open("control_files/"+options.control_name, "w")
    file.write(options_string)
    file.close()


if __name__ == "__main__":
    main()

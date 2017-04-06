#!/bin/bash

# In this file we alter the defaults to generate the CUTEst control files
# for the standard sets of options

# standard methods....
./generate_control_file.py GaussNewton --model 1 \
                                       --maxit 1000 \
                                       --exact_second_derivatives T \
                                       --output_progress_vectors T
./generate_control_file.py Newton -b GaussNewton --model 2
./generate_control_file.py Hybrid -b GaussNewton --model 3

# Newton-Tensor models
./generate_control_file.py NewtonTen2 \
                           --print_level 1 \
                           --maxit 1000 \
                           --model 4 \
                           --type_of_method 2 \
                           --radius_increase 10.0 --radius_reduce 0.1 \
                           --exact_second_derivatives T \
                           --output_progress_vectors T
./generate_control_file.py NewtonTen -b NewtonTen2 --inner_method 1
./generate_control_file.py BaseReg -b NewtonTen2 --inner_method 1
./generate_control_file.py NT_p2_explicit -b BaseReg --inner_method 5
./generate_control_file.py NT_p3_explicit -b BaseReg --inner_method 2
./generate_control_file.py NT_p2_implicit -b BaseReg --inner_method 3
./generate_control_file.py NT_p3_implicit -b BaseReg --inner_method 4

# For the compilation test
./generate_control_file.py TESTSPEC --print_level 3 \
                                    --model 2 \
                                    --maxit 1000 \
                                    --stop_g_absolute 1.0e-7 \
                                    --stop_g_relative 1.0e-7 \
                                    --exact_second_derivatives T \
                                    --output_progress_vectors T

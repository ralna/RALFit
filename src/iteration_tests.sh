#!/bin/bash

# A script to run a sequence of tests
echo "model = 1"
sed -i '5s/.*/1		    ! model/' cutest_control.in
./run_tests.sh results_m1.out

echo "model = 2"
sed -i '5s/.*/2		    ! model/' cutest_control.in
./run_tests.sh results_m2.out

echo "model = 7"
sed -i '5s/.*/7		    ! model/' cutest_control.in
./run_tests.sh results_m7.out

echo "model = 8"
sed -i '5s/.*/8		    ! model/' cutest_control.in
./run_tests.sh results_m8.out

./analyse_iterations.py

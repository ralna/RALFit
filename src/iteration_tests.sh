#!/bin/bash

# A script to run a sequence of tests
echo "model = 1"
sed -i '5s/.*/1                     model used (1=first order, 2=Newton)/' \
    $CUTEST/src/ral_nlls/RAL_NLLS.SPC
./run_tests.sh results_m1.out

echo "model = 2"
sed -i '5s/.*/2                     model used (1=first order, 2=Newton)/' \
    $CUTEST/src/ral_nlls/RAL_NLLS.SPC
./run_tests.sh results_m2.out

echo "model = 7"
sed -i '5s/.*/7                     model used (1=first order, 2=Newton)/' \
    $CUTEST/src/ral_nlls/RAL_NLLS.SPC
./run_tests.sh results_m7.out

echo "model = 8"
sed -i '5s/.*/8                     model used (1=first order, 2=Newton)/' \
    $CUTEST/src/ral_nlls/RAL_NLLS.SPC
./run_tests.sh results_m8.out

echo "gsl"
./run_tests_gsl.sh results_gsl.out

./analyse_iterations.py

#!/bin/bash

# A script to run a sequence of tests
echo "base test (Delta = 100)"
sed -i '7s/.*/0                     relative trust region? (1=relative)/' \
    $CUTEST/src/ral_nlls/RAL_NLLS.SPC
sed -i '9s/.*/100.0                 initial TR radius/' \
    $CUTEST/src/ral_nlls/RAL_NLLS.SPC
./run_tests.sh results_Deltabase.out

# now, let's do the 'real' tests 
sed -i '7s/.*/1                     relative trust region? (1=relative)' \
    $CUTEST/src/ral_nlls/RAL_NLLS.SPC

declare -a lines=( '-3' '0' '3')
for i in "${lines[@]}"
do
    sed -i '9s/.*/1.0E'$i'                TR scaling parameter/' \
    $CUTEST/src/ral_nlls/RAL_NLLS.SPC
    ./run_tests.sh 'results_Delta'$i'.out'
done

#echo "model = 7"
#sed -i '5s/.*/7                     model used (1=first order, 2=Newton)/' \
#    $CUTEST/src/ral_nlls/RAL_NLLS.SPC
#./run_tests.sh results_m7.out#

#echo "model = 8"
#sed -i '5s/.*/8                     model used (1=first order, 2=Newton)/' \
#    $CUTEST/src/ral_nlls/RAL_NLLS.SPC
#./run_tests.sh results_m8.out

#echo "gsl"
#./run_tests_gsl.sh results_gsl.out

./analyse_iterations.py

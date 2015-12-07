#!/bin/bash

# First, make sure the model is Gauss-Newton
sed -i '5s/.*/1                     model used (1=first order, 2=Newton)/' \
    $CUTEST/src/ral_nlls/RAL_NLLS.SPC

# A script to run a sequence of tests
echo "nlls_method = 1"
sed -i '4s/.*/1                     method used (1=dogleg, 2=AINT, 3=More-Sorensen)/' \
    $CUTEST/src/ral_nlls/RAL_NLLS.SPC
./run_tests.sh results_m1.out

echo "nlls_method = 2"
sed -i '4s/.*/2                     method used (1=dogleg, 2=AINT, 3=More-Sorensen)/' \
    $CUTEST/src/ral_nlls/RAL_NLLS.SPC
./run_tests.sh results_m2.out

echo "nlls_method = 3"
sed -i '4s/.*/3                     method used (1=dogleg, 2=AINT, 3=More-Sorensen)/' \
    $CUTEST/src/ral_nlls/RAL_NLLS.SPC
./run_tests.sh results_m3.out

./analyse_subproblem.py

#!/bin/bash

make -C $RAL_NLLS../debug/
RESULT=$?
[ $RESULT -ne 0 ] && exit 1

echo "cutest example..."
cd $RAL_NLLS../test/comparison_tests/cutest/sif/
runcutest --package ral_nlls --architecture pc64.lnx.gfo --decode RAT43.SIF

# run the test suite...
$RAL_NLLS../debug/test/test
RESULT=$?
[ $RESULT -ne 0 ] && exit 1

# run the sample programs...
echo "C example...."
$RAL_NLLS../debug/example/C/nlls_c_example
echo "Fortran example..."
$RAL_NLLS../debug/example/Fortran/nlls_example

#[ $RESULT -ne 0 ] && exit 1
exit 0

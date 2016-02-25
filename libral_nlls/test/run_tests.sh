#!/bin/bash

NLLS_DEBUG=$RAL_NLLS../debug/
NLLS_TEST_SRC=$RAL_NLLS../test/
NLLS_TEST=$NLLS_DEBUG/test/
NLLS_EXAMPLE=$NLLS_DEBUG/example/

make -C $NLLS_DEBUG
RESULT=$?
[ $RESULT -ne 0 ] && exit 1

echo "cutest example..."
cd $NLLS_TEST_SRC/comparison_tests/cutest/sif/
runcutest --package ral_nlls --architecture pc64.lnx.gfo --decode RAT43.SIF | \
    diff $NLLS_TEST_SRC/rat43.output -
RESULT=$?
[ $RESULT -ne 0 ] && exit 2
echo "passed"

# run the test suite...
$NLLS_TEST/test
RESULT=$?
[ $RESULT -ne 0 ] && exit 3

# c tests...
$NLLS_TEST/nlls_c_test
RESULT=$?
[ $RESULT -ne 0 ] && exit 4

# run the sample programs...
echo "C example...."
$NLLS_EXAMPLE/C/nlls_c_example | diff $NLLS_TEST_SRC/nlls_c_example.output -
RESULT=$?
[ $RESULT -ne 0 ] && exit 5
echo "passed"

# $RAL_NLLS../debug/example/C/nlls_c_example
echo "Fortran example..."
$NLLS_EXAMPLE/Fortran/nlls_example | \
    diff <(head -n 164 $NLLS_TEST_SRC/nlls_fortran_example.output) \
    <(head -n 164 -)
# only compare up to the statistics (line 164) -- avoid comparing the printout
# of the timings
RESULT=$?
[ $RESULT -ne 0 ] && exit 6
echo "passed"

#[ $RESULT -ne 0 ] && exit 1
exit 0

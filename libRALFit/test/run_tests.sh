#!/bin/bash

NLLS_BASE=$RAL_NLLS..
NLLS_DEBUG=$NLLS_BASE/build/
#debug/
NLLS_TEST_SRC=$NLLS_BASE/test/
NLLS_TEST=$NLLS_DEBUG/test/
NLLS_EXAMPLE=$NLLS_DEBUG/example/

make clean -C $NLLS_DEBUG

make -C $NLLS_DEBUG
RESULT=$?
[ $RESULT -ne 0 ] && exit 1

cd $NLLS_BASE
rm ral_nlls.so
python setup.py build_ext --inplace
RESULT=$?
[ $RESULT -ne 0 ] && exit 10

echo "cutest example..."
cp $NLLS_TEST_SRC/comparison_tests/control_files/TESTSPEC \
   $NLLS_TEST_SRC/comparison_tests/cutest/sif/RAL_NLLS.SPC
cp $NLLS_TEST_SRC/comparison_tests/cutest/src/ral_nlls/*.f90 \
   $CUTEST/src/ral_nlls/
cd $NLLS_TEST_SRC/comparison_tests/cutest/sif/
runcutest --package ral_nlls --architecture pc64.lnx.gfo --decode RAT43.SIF | \
    diff <(head -n 150 $NLLS_TEST_SRC/rat43.output) <(head -n 150 -)
# skip the final few lines, so as to avoid the timings 
RESULT=$?
[ $RESULT -ne 0 ] && exit 2
echo "passed"

# run the test suite...
echo "fortran tests..."
$NLLS_TEST/nlls_f90_test | diff $NLLS_TEST_SRC/test.output -
RESULT=$?
[ $RESULT -ne 0 ] && exit 3
echo "passed"

# c tests...
echo "C tests..."
$NLLS_TEST/nlls_c_test | diff $NLLS_TEST_SRC/nlls_c_test.output -
RESULT=$?
[ $RESULT -ne 0 ] && exit 4
echo "passed"

# run the sample programs...
echo "C example...."
$NLLS_EXAMPLE/C/nlls_c_example | diff $NLLS_TEST_SRC/nlls_c_example.output -
RESULT=$?
[ $RESULT -ne 0 ] && exit 5
echo "passed"

# $RAL_NLLS../debug/example/C/nlls_c_example
echo "Fortran example..."
$NLLS_EXAMPLE/Fortran/nlls_example | \
    diff $NLLS_TEST_SRC/nlls_fortran_example.output -
RESULT=$?
[ $RESULT -ne 0 ] && exit 6
echo "passed"

echo "Second Fortran example..."
$NLLS_EXAMPLE/Fortran/nlls_example2 | \
    diff $NLLS_TEST_SRC/nlls_fortran_example2.output -
RESULT=$?
[ $RESULT -ne 0 ] && exit 7
echo "passed"

# python!
echo "Python example..."
$NLLS_BASE/example/Python/solve.py | diff $NLLS_TEST_SRC/solve_python.output -
RESULT=$?
[ $RESULT -ne 0 ] && exit 8
echo "passed"

#[ $RESULT -ne 0 ] && exit 1
exit 0

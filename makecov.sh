#!/bin/bash

cd libRALFit
mkdir coverage
cd coverage
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
RESULT=$?
[ $RESULT -ne 0 ] && exit 1
test/nlls_f90_test
RESULT=$?
[ $RESULT -ne 0 ] && exit 2
test/nlls_c_test
RESULT=$?
[ $RESULT -ne 0 ] && exit 3
cd ../src
COV_DIR=../coverage/src/CMakeFiles/ral_nlls.dir
#cp ../coverage/CMakeFiles/ral_nlls.dir/src/ral_nlls_dtrs_double.f90.gcda ral_nlls_dtrs_double.gcda
#cp ../coverage/CMakeFiles/ral_nlls.dir/src/ral_nlls_dtrs_double.f90.gcno ral_nlls_dtrs_double.gcno
#gcov ral_nlls_dtrs_double.f90
#cp ../coverage/CMakeFiles/ral_nlls.dir/src/ral_nlls_symbols.f90.gcno ral_nlls_symbols.gcno
#gcov ral_nlls_symbols.f90
#cp ../coverage/CMakeFiles/ral_nlls.dir/src/ral_nlls_symbols.f90.gcno ral_nlls_double.gcno
#gcov ral_nlls_double.f90
cp ${COV_DIR}/ral_nlls_internal.f90.gcda ral_nlls_internal.gcda
cp ${COV_DIR}/ral_nlls_internal.f90.gcno ral_nlls_internal.gcno
gcov ral_nlls_internal.f90
cp ${COV_DIR}/ral_nlls_workspaces.f90.gcda ral_nlls_workspaces.gcda
cp ${COV_DIR}/ral_nlls_workspaces.f90.gcno ral_nlls_workspaces.gcno
gcov ral_nlls_workspaces.f90
exit 0

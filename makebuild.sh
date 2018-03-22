#!/bin/bash

cd ./libRALFit/
mkdir build
cd build
cmake ..
make 
RESULT=$?
[ $RESULT -ne 0 ] && exit 1
cd test
./nlls_f90_test
RESULT=$?
[ $RESULT -ne 0 ] && exit 2
./nlls_c_test
RESULT=$?
[ $RESULT -ne 0 ] && exit 3
cd ../../

exit 0

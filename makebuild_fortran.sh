#!/bin/bash

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

###########
## build ##
###########

cd $SCRIPTPATH/libRALFit/
mkdir build
cd build
cmake ..
make
make install
RESULT=$?
[ $RESULT -ne 0 ] && exit 1

#######################
## run fortran tests ##
#######################

cd test
./nlls_f90_test
RESULT=$?
[ $RESULT -ne 0 ] && exit 2

exit 0

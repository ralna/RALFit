#!/bin/bash

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

JOBS=${NJOBS:-4}

###########
## build ##
###########

cd $SCRIPTPATH/libRALFit/
mkdir -pv build
cd build

# Make CI/HTML plugin happy
mkdir -pv coverage
echo "For coverage info check GNU gfortran-debug build workspace." > coverage/index.html

echo "Starting build"
echo '$JOBS='"$JOBS"
echo "CWD: `pwd`"

#########################
## Setup a virtualenv  ##
#########################
# TODO REMOVE ONCE THE ALL/PYTHON TARGET IS FIXED
python3 -m venv p3venv
source p3venv/bin/activate

echo "Building configuration: cmake .. ${RALFIT_FLAGS}"
cmake .. ${RALFIT_FLAGS}
make -j${JOBS}
RESULT=$?
[ $RESULT -ne 0 ] && exit 1

make -j${JOBS} install
RESULT=$?
[ $RESULT -ne 0 ] && exit 1

#######################
## run fortran tests ##
#######################

cd test

echo 'Begin Fortran test nlls_f90_test (with RALFIT_UT_CMD_ARGS='$RALFIT_UT_CMD_ARGS')'

# Set traps for  *most* core dumps
trap "exit 31" SIGABRT SIGBUS SIGILL SIGQUIT SIGSEGV 
trap -p

./nlls_f90_test $RALFIT_UT_CMD_ARGS
RESULT=$?
FRESULT=$RESULT
exit $FRESULT
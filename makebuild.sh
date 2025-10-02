#!/bin/bash

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

###########
## build ##
###########

cd $SCRIPTPATH/libRALFit/
mkdir build
cd build

#########################
## Setup a virtualenv  ##
#########################
# TODO REMOVE ONCE THE ALL/PYTHON TARGET IS FIXED
python3 -m venv p3venv
source p3venv/bin/activate

cmake .. ${RALFIT_FLAGS}
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

#################
## run c tests ##
#################
./nlls_c_test &> nlls_c_test.output
RESULT=$?
[ $RESULT -ne 0 ] && exit 3
diff nlls_c_test.output $SCRIPTPATH/libRALFit/test/nlls_c_test.output
RESULT=$?
if [ $RESULT -gt 1 ] 
then
  echo "[C test]: Something is wrong with the diff"
  exit 4
elif [ $RESULT -eq 1 ]
then
  echo "[C test]: Something is wrong with the diff"
  exit 4
else
  echo "** C test passed successfully **"
fi

######################
## run python tests ##
######################

export LD_LIBRARY_PATH=$SCRIPTPATH/libRALFit/build/src/:$LD_LIBRARY_PATH
./nlls_python_test &> nlls_python_test.output
RESULT=$?
if [ $RESULT -ne 0 ]
then
   echo "[Python test]: Failed"
   exit 5
fi
diff nlls_python_test.output $SCRIPTPATH/libRALFit/test/solve_python.output
RESULT=$?
if [ $RESULT -gt 1 ] 
then
  echo "[Python test]: Something is wrong with the diff"
  exit 5
elif [ $RESULT -eq 1 ]
then
  echo "[Python test]: Something is wrong with the diff"
  exit 5
else
  echo "** Python test passed successfully **"
fi

exit 0

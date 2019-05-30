#!/bin/bash

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
cd $SCRIPTPATH/libRALFit/
mkdir build
cd build
cmake ..
make
make install
RESULT=$?
[ $RESULT -ne 0 ] && exit 1
cd test
./nlls_f90_test
RESULT=$?
[ $RESULT -ne 0 ] && exit 2
./nlls_c_test
RESULT=$?
[ $RESULT -ne 0 ] && exit 3
cd $SCRIPTPATH/libRALFit/

######################
## run python tests ##
######################

#cd $SCRIPTPATH/libRALFit/
export LD_LIBRARY_PATH=$SCRIPTPATH/libRALFit/build/src/:$LD_LIBRARY_PATH
#cd test
./nlls_python_test &> nlls_python_test.output
RESULT=$?
[ $RESULT -ne 0 ] && exit 5
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

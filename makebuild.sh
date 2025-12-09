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

echo 'Begin Fortran test nlls_f90_test'

./nlls_f90_test
RESULT=$?
FRESULT=$RESULT
if [ $RESULT -ne 0 ]; then
    echo 'End Fortran test nlls_f90_test: FAILED!'
else
    echo "End Fortran test nlls_f90_test: passed successfully"
fi

# Quick exit if NAGFOR compiler
if [ "$FC" == "nagfor" ]; then
    echo "Note: NAGFOR compiler - skipping C and Python tests"
    exit 0
fi

#################
## run c tests ##
#################

echo 'Begin C test nlls_c_test'

./nlls_c_test > nlls_c_test.output 2> nlls_c_test.stderr
RESULT=$?
CRESULT=$RESULT
if [ $CRESULT -ne 0 ]; then
    echo 'End C test nlls_c_test: FAILED!'
else
    echo "End C test nlls_c_test: passed successfully"
fi

if [ $CRESULT -ne 0 ]; then
    echo -e "\n\nstdout:"
    cat nlls_c_test.output
    echo -e "\n\nstderr:"
    cat nlls_c_test.stderr
fi

diff nlls_c_test.output $SCRIPTPATH/libRALFit/test/nlls_c_test.output
RESULT=$?
if [ $RESULT -ne 0 ]; then
  CRESULT=14 # request exit
  echo "[C test]: Something is wrong with the diff"
  echo "Diff file content:"
  diff -y nlls_c_test.output $SCRIPTPATH/libRALFit/test/nlls_c_test.output
fi

if [ -s nlls_c_test.stderr ]; then
  echo "** C test passed successfully WITH RUNTIME WARNINGS (see nlls_c_test.stderr) **"
  head nlls_c_test.stderr
else
  echo "** C test passed successfully **"
fi

# Stop if either Fortran or C test failed
if [ $FRESULT -ne 0 ] || [ $CRESULT -ne 0 ]; then
    echo "Fortran or C test failed - stopping here."
    exit 2
fi

# Go back
cd $SCRIPTPATH/libRALFit/build

if [[ "$RALFIT_FLAGS" =~ 'CMAKE_BUILD_TYPE=Debug' && "$FC" == "gfortran" ]]; then
        echo "Executing the coverage target"
        make coverage
fi

######################
## run python tests ##
######################

# ASAN may not play nice with python - so disable the python tests
# Infer if debug mode is on
# TODO FIXME re-enable once -DASAN=Off option is added to CMakeLists.txt
if [[ "$RALFIT_FLAGS" =~ 'CMAKE_BUILD_TYPE=Debug' ]]; then
    echo "Note: Debug (ASAN) mode - skipping python tests"
    exit 0
fi

export LD_LIBRARY_PATH=$SCRIPTPATH/libRALFit/build/src/:$LD_LIBRARY_PATH

echo 'Begin Python test nlls_python_test'

../test/nlls_python_test &> nlls_python_test.output
RESULT=$?

echo 'End Python test nlls_python_test: return code' $RESULT

if [ $RESULT -ne 0 ]
then
   echo "[Python test]: Failed"
   cat nlls_python_test.output
   exit 5
fi
diff nlls_python_test.output $SCRIPTPATH/libRALFit/test/solve_python.output 
RESULT=$?
if [ $RESULT -ne 0 ] 
then
  echo "[Python test]: Something is wrong with the diff"
  exit 4
else
  echo "** Python test passed successfully **"
fi

echo "All tests executed. Bye."

exit 0

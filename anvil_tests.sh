#!/bin/bash
# Anvil launch script
# Usage: anvil_tests.sh <compiler> [<extra_flags>]

echo 'AXIS   Building compiler:   $compiler='$compiler
echo 'PARAM  Extra flags:         $CMAKE_EXTRA='$CMAKE_EXTRA

case $compiler in
gfortran|gcc)
   module load gcc/latest
   export CC=gcc
   export F77=gfortran
   export FC=gfortran
   export RALFIT_FLAGS="-DCMAKE_BUILD_TYPE=Release $CMAKE_EXTRA"
   ;;
gfortran-debug|gcc-debug)
   module load gcc/latest
   export CC=gcc
   export F77=gfortran
   export FC=gfortran
   export RALFIT_FLAGS="-DCMAKE_BUILD_TYPE=Debug $CMAKE_EXTRA"
   ;;
ifort|ifx|icx|intel)
   module load intel/latest
   export CC=icx
   export F77=ifx
   export FC=ifx
   export RALFIT_FLAGS="-DCMAKE_BUILD_TYPE=Release $CMAKE_EXTRA"
   ;;
nagfor*) 
   module load gcc/latest
   module load nag/7.2
   export CC=gcc
   export F77=nagfor
   export FC=nagfor
   export RALFIT_FLAGS="-DCMAKE_BUILD_TYPE=Debug $CMAKE_EXTRA"
   ;;
aocc|amd) 
   module load amd
   export CC=clang
   export F77=flang
   export FC=flang
   export RALFIT_FLAGS="-DCMAKE_BUILD_TYPE=Release $CMAKE_EXTRA"
   ;;
aocc-debug|amd-debug) 
   module load amd
   export CC=clang
   export F77=flang
   export FC=flang
   export RALFIT_FLAGS="-DCMAKE_BUILD_TYPE=Debug $CMAKE_EXTRA"
   ;;
*)
    echo "$0 Error: Unknown compiler \"$compiler\"?"
    exit 10
esac


module load cmake/latest
# TODO Selectively use OpenBLAS OR AOCL
case $compiler in
    aocc*)
        export BLAS_LIBRARIES=
        clang --version
        echo '$AOCL_ROOT='$AOCL_ROOT
    ;;
    *)
        module load openblas/latest
        export BLAS_LIBRARIES=-lopenblas
    ;;
esac

echo '$BLAS_LIBRARIES='$BLAS_LIBRARIES

./makebuild.sh

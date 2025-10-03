#module spider
#


echo "Compiler $compiler on target $target"
case $compiler in
gfortran)
   module load gcc/latest
   export CC=gcc
   export F77=gfortran
   export FC=gfortran
   export CFLAGS="-march=native -O3 -Wall -fopenmp"
   export FFLAGS="-march=native -O3 -Wall -fopenmp"
   export RALFIT_FLAGS="-DCMAKE_BUILD_TYPE=Release"
   ;;
gfortran-debug)
   module load gcc/latest
   export CC=gcc
   export F77=gfortran
   export FC=gfortran
   export CFLAGS="-g -O2 -Wall -pedantic -fno-omit-frame-pointer -fopenmp"
   export FFLAGS="-g -O2 -Wall -pedantic -fcheck=all -fbacktrace -fno-omit-frame-pointer -finit-real=nan -finit-integer=-9999 -fopenmp"
   export RALFIT_FLAGS="-DCMAKE_BUILD_TYPE=Debug"
   ;;
ifort|ifx|icx|intel)
   module load intel/latest
   export CC=icx
   export F77=ifx
   export FC=ifx
   export RALFIT_FLAGS="-DCMAKE_BUILD_TYPE=Release"
   ;;
nagfor-debug) 
   module load gcc/latest
   module load nag/7.2
   export CC=gcc
   export F77=nagfor
   export FC=nagfor
   export FFLAGS="-g -nan -C=all -C=undefined -u -ieee=full -kind=unique"
   export RALFIT_FLAGS="-DCMAKE_BUILD_TYPE=Debug"
   ;;
aocc|amd) 
   module load amd
   export CC=clang
   export F77=flang
   export FC=flang
   export CFLAGS="-march=native -O3 -Wall -fopenmp"
   export FFLAGS="-march=native -O3 -Wall -fopenmp"
   export RALFIT_FLAGS="-DCMAKE_BUILD_TYPE=Release"
   ;;
aocc-debug|amd-debug) 
   module load amd
   export CC=clang
   export F77=flang
   export FC=flang
   export CFLAGS="-g -O2 -Wall -pedantic -fno-omit-frame-pointer -fopenmp"
   export FFLAGS="-g -O2 -Wall -pedantic -fno-omit-frame-pointer -fopenmp"
   export RALFIT_FLAGS="-DCMAKE_BUILD_TYPE=Debug"
   ;;
*)
    echo "Unknown compiler $compiler"
    exit 10
esac


module load cmake/latest
# TODO Selectively use openblas OR AOCL
case $compiler in
    aocc*)
        export BLAS_LIBRARIES=
        clang --version
        echo $AOCL_ROOT
    ;;
    *)
        module load openblas/latest
        export BLAS_LIBRARIES=-lopenblas
    ;;
esac


case $compiler in
    nagfor-debug)
	./makebuild_fortran.sh
	;;
    *)
	./makebuild.sh
	;;
esac

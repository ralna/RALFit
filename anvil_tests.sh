#module spider
#
echo "Compiler $compiler on target $target"
case $compiler in
gfortran)
   module load gcc/latest
   export CC=gcc
   export F77=gfortran
   export FC=gfortran
   ;;
gfortran-debug)
   module load gcc/latest
   export CC=gcc
   export F77=gfortran
   export FC=gfortran
   export CFLAGS="-g -O2 -Wall -pedantic -fno-omit-frame-pointer -fopenmp"
   export CXXFLAGS="-g -O2 -Wall -pedantic fno-omit-frame-pointer -fopenmp"
   export FFLAGS="-g -O2 -Wall -pedantic -fcheck=all -fbacktrace -fno-omit-frame-pointer -finit-real=nan -finit-integer=-9999 -fopenmp"
   ;;
ifort)
   module load intel/latest
   export CC=icx
   export F77=ifort
   export FC=ifort
   ;;
nagfor) 
   module load gcc/latest
   module load nag/7.2
   export CC=gcc
   export F77=nagfor
   export FC=nagfor
   ;;
nagfor-debug) 
   module load gcc/latest
   module load nag/7.2
   export CC=gcc
   export F77=nagfor
   export FC=nagfor
   export FFLAGS="-g -nan -C=all -C=undefined -u -ieee=full -kind=unique"
esac


#module load cmake/3.3.1 openblas/0.2.14
module load cmake/latest openblas/latest
export BLAS_LIBRARIES=-lopenblas

case $compiler in
    nagfor-debug)
	./makebuild_fortran.sh
	;;
    *)
	./makebuild.sh
	;;
esac

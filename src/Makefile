#Select compiler and options from below
############################
FC = gfortran
CC = gcc -Wall -g
C99_FLAG = -std=c99

FFLAGS = -Wall -pedantic -fbounds-check -march=native -g -Wunused -fbacktrace -fopenmp -fprofile-arcs -ftest-coverage
#FC = ifort
#FFLAGS = -g -C -static -openmp
#FC = nagfor
#FFLAGS = -dcfuns -C=all -C=undefined -gline -f95

#MKLLOC = /opt/intel/mkl/10.1.0.015/lib/em64t #fox
MKLLOC = /opt/intel/mkl/lib/intel64/         #irons

MKLFLAGS = -lmkl_gf_lp64 -lmkl_sequential -lmkl_core #gfortran, sequential
#MKLFLAGS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core #ifort, sequential
#MKLFLAGS = -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core #gfortran, threaded
#MKLFLAGS = -lmkl_intel_lp64 -lmkl_intel_thread	 -lmkl_core #ifort, threaded

# Get the commands for the libraries:
LIBS = -L$(MKLLOC) $(MKLFLAGS) 

# Compile f77 code with f90 compiler
F77 = $(FC)

# What to build
all: driver

DLAPACK = $(LAPACKDIR)/dgels.f

#DDEPS = 

#COMMON90 =

#SDEPS90 = 

#DDEPS90 = 


driver.o :  example_module.o nlls_module.o
cdriver.o : example_module.o nlls_module.o
#nlls_module.o : common90.o ddeps90.o ddeps.o  

example_module.o : nlls_module.o

driver: driver.o example_module.o nlls_module.o  # common90.o ddeps90.o ddeps.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

cdriver: cdriver.o example_module.o ral_nlls_ciface.o nlls_module.o # common90.o ddeps90.o ddeps.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

BLASDIR=$(TOPDIR)/blas
LAPACKDIR=$(TOPDIR)/lapack

# These two rules tell us how to compile module files from source files
%.o: %.c
	$(CC) $(CFLAGS) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<
%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

# The following rules build a source file from a variable
common.f: $(COMMON)  $(COMMONBLAS) $(COMMONLAPACK)
	cat $^ > $@
common90.f90: $(COMMON90)
	cat $^ > $@
sdeps.f: $(SDEPS) $(SBLAS) $(SLAPACK)
	cat $^ > $@
ddeps.f: $(DDEPS) $(DBLAS) $(DLAPACK)
	cat $^ > $@
sdeps90.f90: $(SDEPS90)
	cat $^ > $@
ddeps90.f90: $(DDEPS90)
	cat $^ > $@

clean:
	-rm *.o *.mod
	-rm *.gcno *.gcna *.gcda
	-rm metis.f *deps.f *deps90.f90 common.f common90.f90

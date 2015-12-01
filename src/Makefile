############################
#Select compiler and options from below
############################
FC = gfortran
CC = gcc -Wall -g
C99_FLAG = -std=c99
FFLAGS = -Wall -pedantic -fbounds-check -march=native -g -Wunused -fbacktrace -fopenmp -fprofile-arcs -ftest-coverage
############################
#FC = ifort
#FFLAGS = -g -C -static -openmp
############################
#FC = nagfor
#FFLAGS = -dcfuns -C=all -C=undefined -gline -f95
############################

# BLAS
MKL_DIR = /opt/intel/mkl/lib/intel64/
LIBS = -L$(MKL_DIR) -lmkl_gf_lp64 -lmkl_sequential -lmkl_core

# CUTEST
CUTESTLOC = $(CUTEST)/objects/pc64.lnx.gfo/double/
LIBS += -L$(CUTESTLOC) -lcutest	

# GSL
GSL_LIBS = -lgslcblas -lgsl -lm

# What to build
all: cdriver driver libral_nlls.a test

driver: driver.o example_module.o libral_nlls.a
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

test: nlls_test.o example_module.o libral_nlls.a
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

#CUTEST_EXTRA_OBJ = ELFUN.o GROUP.o RANGE.o
#cutest: run_cutest.o libral_nlls.a $(CUTEST_EXTRA_OBJ)
#	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

libral_nlls.a: nlls_module.o ral_nlls_ciface.o
	$(AR) rcs $@ $^

gen_problem:
	sifdecoder $(SIFFILE).SIF

cdriver: cdriver.o example_module.o libral_nlls.a
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS) $(CLIBS)

lm_example: lm_example.o libral_nlls.a
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS) $(CLIBS)

# Tell the makefile where to look for sif files
SIFFILESLOC = ../cutest/sif
SIFFILE = ARGAUSS


# These two rules tell us how to compile module files from source files
%.o: %.c
	$(CC) $(CFLAGS) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<
%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	-rm libral_nlla.a driver cdriver lm_example
	-rm *.o *.mod
	-rm *.gcno *.gcna *.gcda

driver.o :  example_module.o nlls_module.o
test.o :  example_module.o nlls_module.o
cdriver.o : example_module.o nlls_module.o
lm_example.o : nlls_module.o
run_cutest.o: nlls_module.o
example_module.o : nlls_module.o

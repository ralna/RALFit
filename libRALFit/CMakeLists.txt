# CMake project file for libral_nlls
cmake_minimum_required (VERSION 2.6)
enable_language (C Fortran)

project (RALFit)

# Option to compile our own LAPACK
option(CompileMiniLAPACK "Compile our own LAPACK (discouraged)" Off)

option(TestCoverage "Perform code coverage analysis" Off)

option(UseMKL "Use Intel's MKL library" Off)

# Set the location for cmake module files
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

set(CMAKE_EXE_LINKER_FLAGS "-Wl,-export-dynamic")

# Set the source files that we're going to compile
set ( SRC_FILES
        src/ral_nlls_double.f90
	src/ral_nlls_internal.f90
	src/ral_nlls_dtrs_double.f90
	src/ral_nlls_symbols.f90
	src/ral_nlls_ciface.f90
	src/ral_nlls_workspaces.f90
	)	

if (UseMKL)
   find_package(MKL REQUIRED)
   if ( NOT MKL_FOUND )
      message( FATAL_ERROR "MKL required, but not found! Check $MKLROOT or $INTEL environment variables are correctly set.")
  endif ( NOT MKL_FOUND ) 				
  message( STATUS "The MKL library was found at" ${MKL_LIBRARIES} )
  set(LIBS ${LIBS} ${MKL_LIBRARIES}) # ${MKL_INTERFACE_LIBRARY} ${MKL_SEQUENTIAL_LAYER_LIBRARY} ${MKL_CORE_LIBRARY})	
  include_directories(${MKL_INCLUDE_DIRS})
  set(CompileMiniLAPACK Off)
endif(UseMKL)

# see if the build type is debug.  First, convert to lower
STRING(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_LOWER)
# now do the test
if(CMAKE_BUILD_TYPE_LOWER MATCHES debug)
  message("WARNING: Debug version - not to be used in production code.")
  # compile lapack from source...
  set(CompileMiniLAPACK On)
  # set the test coverage switch
  set(TestCoverage On)
elseif(CMAKE_BUILD_TYPE_LOWER MATCHES profile)
      message("WARNING: Debug version - not to be used in production code.")	
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -pg -g -O3")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg -g -O3")
  set(CMAKE_LDFLAGS "${CMAKE_Fortran_FLAGS} -pg -g -O3")		    
else()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
endif()


# Set flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -fPIC")
#set(CMAKE_Fortran_COMPILER nagfor)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")

set( CMAKE_VERBOSE_MAKEFILE on )		
if(TestCoverage)
   set( CMAKE_VERBOSE_MAKEFILE on )
   INCLUDE(CodeCoverage)
   set(CMAKE_C_FLAGS_DEBUG "-pg -g -O0 -Wall -pedantic -fprofile-arcs -ftest-coverage")
   set(CMAKE_Fortran_FLAGS_DEBUG "-pg -g -O0 -Wall -pedantic -fcheck=all -fbacktrace -finit-real=nan -finit-integer=-9999 -fprofile-arcs -ftest-coverage")
   set(CMAKE_LDFLAGS_DEBUG "-g -fprofile-arcs -ftest-coverage")
   SETUP_TARGET_FOR_COVERAGE(test_coverage test/test coverage)
endif(TestCoverage)


if((NOT CompileMiniLAPACK) AND (NOT UseBLAS))
  # try to find lapack..
  find_package(LAPACK)
  if(NOT LAPACK_FOUND)
    # couldn't find lapack...compile locally
    set(CompileMiniLAPACK On)
  endif()
endif()

# LAPACK routines
if(CompileMiniLAPACK)
   message("WARNING: Compiling our own LAPACK library.")
   message("WARNING: This should be avoided as it is likely to be less optimised than a")
   message("WARNING: vendor supplied version, and may also be out of date.")
   add_library(ral_nlls_lapack
         lapack/dbdsqr.f lapack/dgebd2.f lapack/dgebrd.f lapack/dgebrd.f
         lapack/dgebrd.f lapack/dgelq2.f lapack/dgelqf.f lapack/dgels.f
         lapack/dgeqr2.f lapack/dgeqrf.f lapack/dgesvd.f lapack/dgesv.f
         lapack/dgetrf2.f lapack/dgetrf.f lapack/dgetrs.f lapack/dggbak.f
         lapack/dggbal.f lapack/dggev.f lapack/dggev.f lapack/dgghrd.f
         lapack/dhgeqz.f lapack/disnan.f lapack/dlabad.f lapack/dlabrd.f
         lapack/dlacpy.f lapack/dladiv.f lapack/dlae2.f lapack/dlaebz.f
         lapack/dlaev2.f lapack/dlag2.f lapack/dlagtf.f lapack/dlagts.f
         lapack/dlaisnan.f lapack/dlaln2.f lapack/dlange.f lapack/dlanhs.f
         lapack/dlanst.f lapack/dlansy.f lapack/dlapy2.f lapack/dlapy3.f
         lapack/dlarfb.f lapack/dlarf.f lapack/dlarfg.f lapack/dlarft.f
         lapack/dlarnv.f lapack/dlartg.f lapack/dlaruv.f lapack/dlas2.f
         lapack/dlascl.f lapack/dlaset.f lapack/dlasq1.f lapack/dlasq2.f
         lapack/dlasq3.f lapack/dlasq4.f lapack/dlasq5.f lapack/dlasq6.f
         lapack/dlasr.f lapack/dlasrt.f lapack/dlassq.f lapack/dlasv2.f
         lapack/dlaswp.f lapack/dlatrd.f lapack/dorg2l.f lapack/dorg2r.f
         lapack/dorgbr.f lapack/dorgl2.f lapack/dorglq.f lapack/dorgql.f
         lapack/dorgqr.f lapack/dorgtr.f lapack/dorm2l.f lapack/dorm2r.f
         lapack/dormbr.f lapack/dorml2.f lapack/dormlq.f lapack/dormql.f 
         lapack/dormqr.f lapack/dormtr.f lapack/dposv.f lapack/dpotrf2.f
         lapack/dpotrf.f lapack/dpotrs.f lapack/dstebz.f lapack/dstein.f
         lapack/dsteqr.f lapack/dsterf.f lapack/dsyev.f lapack/dsyevx.f
         lapack/dsytd2.f lapack/dsytrd.f lapack/dtgevc.f lapack/dtrtrs.f
         lapack/ieeeck.f lapack/iladlc.f lapack/iladlr.f lapack/ilaenv.f
         lapack/iparmq.f lapack/dlamch.f)
   set(LAPACK_LIBRARIES ral_nlls_lapack)
   find_package(BLAS REQUIRED)
endif()

message( STATUS BLAS found: ${BLAS_LIBRARIES} )

# Main libral_nlls library
include_directories(include ${CMAKE_BINARY_DIR}) # BINARY_DIR for generated .mod
add_library (ral_nlls SHARED ${SRC_FILES} ) 
target_link_libraries(ral_nlls ${LAPACK_LIBRARIES} -lm -ldl)


# Tests and examples
add_subdirectory (test)
add_subdirectory (example/C)
add_subdirectory (example/Fortran)
#add_subdirectory (example/C++)

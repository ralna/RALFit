# Copyright (c) 2016, The Science and Technology Facilities Council (STFC)
# All rights reserved.
# Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
# CMake project file for libral_nlls
cmake_minimum_required (VERSION 3.1)

project (RALFit LANGUAGES Fortran C)

# Request to run pre-processor as well
if(NOT WIN32)
  set(CMAKE_Fortran_FLAGS -cpp)
else()
  set(CMAKE_Fortran_PREPROCESS ON)
endif()

# Option to compile our own LAPACK
option(CompileMiniLAPACK "Compile our own LAPACK (discouraged)" Off)
option(TestCoverage "Perform code coverage analysis" Off)

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
        src/ral_nlls_printing.f90
        src/ral_nlls_types.F90
        src/ral_nlls_fd.f90
	)

# Support for CMAKE_BUILD_TYPE={Debug,Release,Profile}
if(CMAKE_BUILD_TYPE MATCHES Debug)
  message("WARNING: Debug version - not to be used in production code.")
  # compile lapack from source...
  set(CompileMiniLAPACK On)
  # set the test coverage switch
  set(TestCoverage On)
elseif(CMAKE_BUILD_TYPE MATCHES Profile)
  message("WARNING: Debug version - not to be used in production code.")
endif()


# Set flags
set(CMAKE_POSITION_INDEPENDENT_CODE ON) # SET fPIC
set(CMAKE_C_STANDARD 99) # SET -std=c99

if(TestCoverage)
   set( CMAKE_VERBOSE_MAKEFILE on )
   INCLUDE(CodeCoverage)
   set(CMAKE_C_FLAGS_DEBUG "-pg -g -O0 -Wall -pedantic -fprofile-arcs -ftest-coverage -fsanitize=address -fsanitize=leak")
   set(CMAKE_Fortran_FLAGS_DEBUG "-pg -g -O0 -Wall -pedantic -fcheck=all -fbacktrace -finit-real=nan -finit-integer=-9999 -fprofile-arcs -ftest-coverage -fsanitize=address -fsanitize=leak")
   set(CMAKE_LDFLAGS_DEBUG "-g -fprofile-arcs -ftest-coverage -fsanitize=address -fsanitize=leak")
   SETUP_TARGET_FOR_COVERAGE(test_coverage test/nlls_f90_test coverage)
endif(TestCoverage)

# BLAS and LAPACK

if(NOT CompileMiniLAPACK)
  # try to find lapack..
  find_package(LAPACK)
  if(NOT LAPACK_FOUND)
    # couldn't find lapack...compile locally
    set(CompileMiniLAPACK On)
  else()
    message( STATUS "The LAPACK library was found at" ${LAPACK_LIBRARIES} )
    set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
  endif()
endif()

# LAPACK routines
if(CompileMiniLAPACK)
   message("*****************************WARNING**********************************")
   message("* Compiling our own LAPACK library.                                  *")
   message("* This should be avoided as it is likely to be less optimised than a *")
   message("* vendor supplied version, and may also be out of date.              *")
   message("**********************************************************************")
   # note that blas is a dependency of lapack, so blas only explicitly needed here,
   # if we've built our own lapack
   find_package(BLAS REQUIRED)
   if ( NOT BLAS_FOUND )
     message( FATAL_ERROR "BLAS required, but not found!")
   endif ( NOT BLAS_FOUND )
   message( STATUS "The BLAS library was found at" ${BLAS_LIBRARIES} )
   set(LIBS ${LIBS} ${BLAS_LIBRARIES})
   add_subdirectory(lapack)
   set(LAPACK_LIBRARIES ral_nlls_lapack)
   set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
endif()

# Python interface
find_package (PythonInterp)

if (PYTHONINTERP_FOUND)
   message( STATUS "looking for file in ${CMAKE_CURRENT_SOURCE_DIR}")
   set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
   set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
   set(OUTPUT   "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")

   configure_file(${SETUP_PY_IN} ${SETUP_PY})

   execute_process (
     COMMAND which ${PYTHON_EXECUTABLE}
     OUTPUT_VARIABLE outVar
     )
   message( STATUS ${outVar} )

   add_custom_command(OUTPUT ${OUTPUT}
   		      COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} build
		      COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT})


   message( STATUS "python found -- writing to ${CMAKE_CURRENT_BINARY_DIR}")

   add_custom_target(target ALL DEPENDS ${OUTPUT} ral_nlls)

   install(CODE "execute_process(
   	    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --user
	    )"
	   )
endif()

# Main libral_nlls library
add_subdirectory(src)


# Tests and examples
add_subdirectory (test)
add_subdirectory (example/C)
add_subdirectory (example/Fortran)
add_subdirectory (example/Python)
#add_subdirectory (example/C++)

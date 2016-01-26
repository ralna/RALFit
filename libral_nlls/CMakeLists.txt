# CMake project file for libral_nlls
cmake_minimum_required (VERSION 2.6)
enable_language (C Fortran)

project (libral_nlls)


# Set the location for cmake module files
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# Set flags
set(CMAKE_C_FLAGS "-g -O2 -std=c99")
INCLUDE(CodeCoverage)
set(CMAKE_C_FLAGS_DEBUG "-g -O0 -fprofile-arcs -ftest-coverage")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fprofile-arcs -ftest-coverage")
set(CMAKE_LDFLAGS_DEBUG "-fprofile-arcs -ftest-coverage")
SETUP_TARGET_FOR_COVERAGE(test_coverage test/test coverage)
#
#set(CMAKE_C_FLAGS_DEBUG "-g -O0 -std=c99 -Wall -W -fprofile-arcs -ftest-coverage")
#set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -Wall -W -fprofile-arcs -ftest-coverage")
#set(CMAKE_LDFLAGS_DEBUG "-fprofile-arcs -ftest-coverage")

# Main libral_nlls library
include_directories(include ${CMAKE_BINARY_DIR}) # BINARY_DIR for generated .mod
add_library (ral_nlls src/ral_nlls_double.f90 src/ral_nlls_internal.f90 src/ral_nlls_dtrs_double.f90 src/ral_nlls_symbols.f90 src/ral_nlls_ciface.f90)

# Tests and examples
find_package(LAPACK REQUIRED)
add_subdirectory (test)
add_subdirectory (example/C)
add_subdirectory (example/Fortran)

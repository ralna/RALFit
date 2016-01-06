# CMake project file for libral_nlls
cmake_minimum_required (VERSION 2.6)
enable_language (C Fortran)
project (libral_nlls)

# Main libral_nlls library
include_directories(include ${CMAKE_BINARY_DIR}) # BINARY_DIR for generated .mod
add_library (ral_nlls src/nlls_module.f90 src/ral_nlls_ciface.f90)

# Tests and examples
find_package(LAPACK REQUIRED)
add_subdirectory (test)
add_subdirectory (example)

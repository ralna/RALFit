# CMake project file for libral_nlls
cmake_minimum_required (VERSION 2.6)
enable_language (Fortran)
project (ral_nlls)

# Main libral_nlls library
add_library (libral_nlls src/nlls_module.f90 src/ral_nlls_ciface.f90)

# Test code
add_executable (test test/nlls_test.f90)
target_link_libraries(test libral_nlls)

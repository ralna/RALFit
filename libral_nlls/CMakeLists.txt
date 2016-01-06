# CMake project file for libral_nlls
cmake_minimum_required (VERSION 2.6)
enable_language (Fortran)
project (libral_nlls)

# Main libral_nlls library
add_library (ral_nlls src/nlls_module.f90 src/ral_nlls_ciface.f90)

# Test code
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
add_executable (test test/nlls_test.f90 test/example_module.f90)
target_link_libraries(test ral_nlls ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# build instructions for test codes

add_executable (test nlls_test.f90 example_module.f90)
target_link_libraries(test ral_nlls ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_executable (cdriver cdriver.c)
target_link_libraries(cdriver ral_nlls ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

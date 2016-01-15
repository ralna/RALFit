# Examples build

add_executable (lm_example lm_example.c expfit.c)
target_link_libraries(lm_example ral_nlls ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

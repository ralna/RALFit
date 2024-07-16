# Copyright (c) 2019, The Science and Technology Facilities Council (STFC)
# All rights reserved.
# Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
# Examples build

# lm_example requires GSL libraries!
#add_executable (lm_example lm_example.c expfit.c)
#target_link_libraries(lm_example ral_nlls ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_executable (nlls_c_example nlls_example.c)
target_link_libraries(nlls_c_example ral_nlls -lm)

add_executable (nlls_c_example2 nlls_example2.c)
target_link_libraries(nlls_c_example2 ral_nlls -lm)

add_executable (nlls_c_example_box nlls_example_box.c)
target_link_libraries(nlls_c_example_box ral_nlls -lm)

add_executable (nlls_c_example_fd nlls_example_fd.c)
target_link_libraries(nlls_c_example_fd ral_nlls -lm)

add_executable (nlls_c_example_chk nlls_example_chk.c)
target_link_libraries(nlls_c_example_chk ral_nlls -lm)

add_executable (tracks tracks.c)
target_link_libraries(tracks ral_nlls -lm)
# Examples build

add_executable (nlls_example nlls_example.f90)
target_link_libraries(nlls_example ral_nlls)

add_executable (nlls_example2 nlls_example2.f90)
target_link_libraries(nlls_example2 ral_nlls)

add_executable (lanczos Lanczos.f90)
target_link_libraries(lanczos ral_nlls)

add_executable (LanczosBox LanczosBox.f90)
target_link_libraries(LanczosBox ral_nlls)

add_executable (stp_x0sol stp_x0sol.f90)
target_link_libraries(stp_x0sol ral_nlls)

#add_executable ( nlls_example_iterate nlls_example_iterate.f90)
#target_link_libraries(nlls_example_iterate ral_nlls -lm)

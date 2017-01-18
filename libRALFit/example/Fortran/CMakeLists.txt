# Examples build

add_executable (nlls_example nlls_example.f90)
target_link_libraries(nlls_example ral_nlls)

add_executable (nlls_example2 nlls_example2.f90)
target_link_libraries(nlls_example2 ral_nlls)

#add_executable ( nlls_example_iterate nlls_example_iterate.f90)
#target_link_libraries(nlls_example_iterate ral_nlls -lm)

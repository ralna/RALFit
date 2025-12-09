# Copyright (c) 2019, The Science and Technology Facilities Council (STFC)
# All rights reserved.
# build instructions for test codes

add_executable (nlls_f90_test unit_test_mod.F90 nlls_test.F90 )
target_link_libraries(nlls_f90_test ral_nlls -lm)

add_executable (cdriver cdriver.c)
target_link_libraries(cdriver ral_nlls -lm)

add_executable (nlls_c_test nlls_c_test.c)
target_link_libraries(nlls_c_test ral_nlls -lm)


if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
set_source_files_properties(unit_test_mod.F90 PROPERTIES COMPILE_FLAGS "$<$<CONFIG:DEBUG>:-Wno-unused-dummy-argument>")
# elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
set_source_files_properties(unit_test_mod.F90 PROPERTIES COMPILE_FLAGS "$<$<CONFIG:DEBUG>:-w=uda>")
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/nlls_python_test
	${CMAKE_CURRENT_BINARY_DIR}/nlls_python_test
 	COPYONLY)

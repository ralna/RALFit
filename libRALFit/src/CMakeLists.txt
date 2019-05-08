add_library (ral_nlls SHARED
	ral_nlls_double.f90
	ral_nlls_internal.f90
	ral_nlls_dtrs_double.f90
	ral_nlls_symbols.f90
	ral_nlls_ciface.f90
	ral_nlls_workspaces.f90
	ral_nlls_printing.f90
	) 

target_link_libraries(ral_nlls ${LIBS} m ${CMAKE_DL_LIBS})
target_include_directories(ral_nlls PUBLIC ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/src)

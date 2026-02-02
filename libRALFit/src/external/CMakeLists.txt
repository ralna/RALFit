if(SINGLE_PRECISION)
    set(precision "single")
else()
    set(precision "double")
endif()
add_library(lsqr_reverse_${precision} OBJECT lsqr.f90)
add_library(dct_module_${precision} OBJECT dct.f90)

target_link_libraries(dct_module_${precision} fftw3)

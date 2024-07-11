# Copyright (c) 2016, The Science and Technology Facilities Council (STFC)
# All rights reserved.
# Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.

set( source_files
  ral_nlls_double.f90
  ral_nlls_internal.f90
  ral_nlls_dtrs_double.f90
  ral_nlls_symbols.f90
  ral_nlls_ciface.f90
  ral_nlls_workspaces.f90
  ral_nlls_printing.f90
  ral_nlls_bounds.f90
  nag_export_mod.f90
  ral_nlls_types.F90
  ral_nlls_fd.f90
  )

add_library (ral_nlls-static STATIC ${source_files})
add_library (ral_nlls SHARED ${source_files})

target_link_libraries(ral_nlls ${LIBS} m ${CMAKE_DL_LIBS})
target_include_directories(ral_nlls PUBLIC ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/src)

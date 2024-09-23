# Copyright (c) 2016, The Science and Technology Facilities Council (STFC)
# All rights reserved.
# Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.

set( source_files
  ral_nlls_types.F90
  ral_nlls_symbols.F90
  ral_nlls_dtrs.F90
  ral_nlls_workspaces.F90
  ral_nlls_fd.F90
  nag_export_mod.F90
  ral_nlls_printing.F90
  ral_nlls_bounds.F90
  ral_nlls_internal.F90
  ral_nlls.F90
  ral_nlls_ciface.F90
  CACHE INTERNAL ""
  )

# compiler warning suppression list
set_source_files_properties(ral_nlls_ciface.F90 ral_nlls_fd.F90 ral_nlls_internal.F90 ral_nlls_workspaces.F90 PROPERTIES COMPILE_FLAGS "$<$<CONFIG:DEBUG>:-Wno-unused-dummy-argument>")
if (SUPPRESS_INTERNAL_FIXMES)
message(AUTHOR_WARNING "Suppressing known compiler warnings for some source files")
# False positives
set_source_files_properties(ral_nlls_internal.F90 PROPERTIES COMPILE_FLAGS "$<$<CONFIG:DEBUG>:-Wno-unused-variable -Wno-unused-dummy-argument>")
# Warning: ‘dtrs_theta_derivs’ defined but not used [-Wunused-function]
set_source_files_properties(ral_nlls_dtrs.F90 PROPERTIES COMPILE_FLAGS "$<$<CONFIG:DEBUG>:-Wno-unused-function>")
endif(SUPPRESS_INTERNAL_FIXMES)

################################################################################
# Build manually dependecies 
set(OBJECT_LIST "")
 foreach(file ${source_files})
    string(REPLACE ".F90" "_lib" obj_name ${file})
    add_library(${obj_name} OBJECT ${file})
    list(APPEND OBJECT_LIST ${obj_name})
 endforeach()
  
 add_dependencies(ral_nlls_dtrs_lib ral_nlls_types_lib ral_nlls_symbols_lib)
 add_dependencies(ral_nlls_workspaces_lib ral_nlls_dtrs_lib)
 add_dependencies(ral_nlls_printing_lib ral_nlls_workspaces_lib)
 add_dependencies(ral_nlls_fd_lib ral_nlls_workspaces_lib ral_nlls_printing_lib)
 add_dependencies(nag_export_mod_lib ral_nlls_workspaces_lib)
 add_dependencies(ral_nlls_bounds_lib ral_nlls_workspaces_lib)

 add_dependencies(ral_nlls_internal_lib ral_nlls_fd_lib nag_export_mod_lib ral_nlls_bounds_lib)

 add_dependencies(ral_nlls_lib ral_nlls_internal_lib)
 add_dependencies(ral_nlls_ciface_lib ral_nlls_lib ral_nlls_internal_lib)
################################################################################

add_library (ral_nlls-static STATIC ${source_files})
add_library (ral_nlls SHARED ${source_files})

add_dependencies(ral_nlls-static ral_nlls_lib ral_nlls_ciface_lib)
add_dependencies(ral_nlls ral_nlls-static)

target_link_libraries(ral_nlls ${LIBS} m ${CMAKE_DL_LIBS})
target_include_directories(ral_nlls PUBLIC ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/src)

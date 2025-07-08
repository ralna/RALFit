# Copyright (c) 2016, The Science and Technology Facilities Council (STFC)
# All rights reserved.
# Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.

# Adapted to build dual precision objects for AOCL-DA

set(source_files
    ral_nlls_types.F90
    ral_nlls_symbols.F90
    ral_nlls_dtrs.F90
    ral_nlls_workspaces.F90
    ral_nlls_printing.F90
    ral_nlls_bounds.F90
    nag_export_mod.F90
    ral_nlls_fd.F90
    ral_nlls_internal.F90
    ral_nlls.F90
    ral_nlls_ciface.F90
    CACHE INTERNAL "source_files" FORCE)

# Suppress warnings ############
set_source_files_properties(ral_nlls_ciface.F90 ral_nlls_fd.F90 ral_nlls_internal.F90 ral_nlls_workspaces.F90 PROPERTIES COMPILE_FLAGS "$<$<CONFIG:DEBUG>:-Wno-unused-dummy-argument>")
if (SUPPRESS_INTERNAL_FIXMES)
message(AUTHOR_WARNING "Suppressing known compiler warnings for some source files")
# False positives
set_source_files_properties(ral_nlls_internal.F90 PROPERTIES COMPILE_FLAGS "$<$<CONFIG:DEBUG>:-Wno-unused-variable -Wno-unused-dummy-argument>")
endif(SUPPRESS_INTERNAL_FIXMES)
# End Supress warnings #########

# Generate target name based on filename (remove file extension and convert to
# uppercase)
function(gen_name filename precision result)
  get_filename_component(name_without_extension ${filename} NAME_WE)
  string(TOUPPER ${name_without_extension} uppercase_name)
  set(${result}
      "${uppercase_name}_${precision}"
      PARENT_SCOPE)
endfunction()

function(generate_ral_nlls_modules precision obj_list)
  set(in_obj_list "")
  foreach(file ${source_files})
    gen_name(${file} ${precision} target)

    add_library(${target} OBJECT ${file})
    if(WIN32 AND CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
      include_directories(" -I${CMAKE_CURRENT_SOURCE_DIR} ")
    endif()
    if(${precision} STREQUAL "single")
      set(local_defs "-DSINGLE_PRECISION")
      target_compile_definitions(${target} PRIVATE ${local_defs})
    endif()
    list(APPEND in_obj_list $<TARGET_OBJECTS:${target}>)
  endforeach()
  set(obj_list
      ${in_obj_list}
      PARENT_SCOPE)

  # Build explicit dependency tree
  add_dependencies("RAL_NLLS_DTRS_${precision}" "RAL_NLLS_TYPES_${precision}"
                   "RAL_NLLS_SYMBOLS_${precision}")
  add_dependencies("RAL_NLLS_WORKSPACES_${precision}"
                   "RAL_NLLS_DTRS_${precision}")
  add_dependencies("NAG_EXPORT_MOD_${precision}"
                   "RAL_NLLS_WORKSPACES_${precision}")
  add_dependencies("RAL_NLLS_PRINTING_${precision}"
                   "RAL_NLLS_WORKSPACES_${precision}")
  add_dependencies("RAL_NLLS_BOUNDS_${precision}"
                   "RAL_NLLS_WORKSPACES_${precision}")
  add_dependencies(
    "RAL_NLLS_FD_${precision}" "RAL_NLLS_WORKSPACES_${precision}"
    "RAL_NLLS_PRINTING_${precision}")
  add_dependencies(
    "RAL_NLLS_INTERNAL_${precision}" "RAL_NLLS_PRINTING_${precision}"
    "RAL_NLLS_BOUNDS_${precision}" "RAL_NLLS_FD_${precision}"
    "NAG_EXPORT_MOD_${precision}" "RAL_NLLS_WORKSPACES_${precision}")
  add_dependencies("RAL_NLLS_${precision}" "RAL_NLLS_INTERNAL_${precision}")
  add_dependencies("RAL_NLLS_CIFACE_${precision}" "RAL_NLLS_${precision}")

  # Replicate tree to be compatible with Ninja
  target_link_libraries("RAL_NLLS_DTRS_${precision}" "RAL_NLLS_TYPES_${precision}"
                        "RAL_NLLS_SYMBOLS_${precision}")
  target_link_libraries("RAL_NLLS_WORKSPACES_${precision}"
                        "RAL_NLLS_DTRS_${precision}")
  target_link_libraries("NAG_EXPORT_MOD_${precision}"
                        "RAL_NLLS_WORKSPACES_${precision}")
  target_link_libraries("RAL_NLLS_PRINTING_${precision}"
                        "RAL_NLLS_WORKSPACES_${precision}")
  target_link_libraries("RAL_NLLS_BOUNDS_${precision}"
                        "RAL_NLLS_WORKSPACES_${precision}")
  target_link_libraries(
    "RAL_NLLS_FD_${precision}" "RAL_NLLS_WORKSPACES_${precision}"
    "RAL_NLLS_PRINTING_${precision}")
  target_link_libraries(
    "RAL_NLLS_INTERNAL_${precision}" "RAL_NLLS_PRINTING_${precision}"
    "RAL_NLLS_BOUNDS_${precision}" "RAL_NLLS_FD_${precision}"
    "NAG_EXPORT_MOD_${precision}" "RAL_NLLS_WORKSPACES_${precision}")
  target_link_libraries("RAL_NLLS_${precision}" "RAL_NLLS_INTERNAL_${precision}")
  target_link_libraries("RAL_NLLS_CIFACE_${precision}" "RAL_NLLS_${precision}")
endfunction()

# Post ################
set(nlls_obj_list "")
if(SINGLE_PRECISION)
    # Generate single-precision library
    generate_ral_nlls_modules("single" obj_list)
else()
    # Generate double-precision library
    generate_ral_nlls_modules("double" obj_list)
endif()
add_library (ral_nlls-static STATIC ${obj_list})
add_library (ral_nlls SHARED ${obj_list})
target_link_libraries(ral_nlls ${LIBS} m ${CMAKE_DL_LIBS})
target_include_directories(ral_nlls PUBLIC ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/src)
# End Post ################

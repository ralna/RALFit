# Copyright (c) 2016, The Science and Technology Facilities Council (STFC)
# All rights reserved.
# Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
# CMake project file for libral_nlls
cmake_minimum_required (VERSION 3.18)

project (RALFit VERSION 2.0 LANGUAGES Fortran C)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
endif()

# Request to run pre-processor
set(CMAKE_Fortran_PREPROCESS ON)

# Option to compile our own LAPACK
option(CompileMiniLAPACK "Compile our own LAPACK (discouraged)" Off)
option(TestCoverage "Perform code coverage analysis" Off)

# Export Release Information automatically
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set the location for cmake module files
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

set(CMAKE_EXE_LINKER_FLAGS "-Wl,-export-dynamic")

# Set the source files that we're going to compile
set ( SRC_FILES
       src/ral_nlls_double.f90
       src/ral_nlls_internal.f90
       src/ral_nlls_dtrs_double.f90
       src/ral_nlls_symbols.f90
       src/ral_nlls_ciface.f90
       src/ral_nlls_workspaces.f90
       src/ral_nlls_printing.f90
       src/ral_nlls_types.F90
       src/ral_nlls_fd.f90
	)

if(CMAKE_BUILD_TYPE MATCHES Release)
  set(TestCoverage Off)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_C_FLAGS_RELEASE "-O2 -Wall -pedantic -march=native -fno-fast-math -fopenmp -fno-omit-frame-pointer")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -Wall -pedantic -march=native -fopenmp -fno-omit-frame-pointer")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
    set(CMAKE_C_FLAGS_RELEASE "-O2 -Wall -pedantic -march=native -fno-fast-math -fopenmp -fno-omit-frame-pointer")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -Wall -pedantic -march=native -fopenmp -fno-omit-frame-pointer")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
    set(CMAKE_C_FLAGS_RELEASE "-O2")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  else()
    set(CMAKE_C_FLAGS_RELEASE "-O2")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  endif()
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
  message("WARNING: Debug version - not to be used in production code.")
  # compile lapack from source...
  set(CompileMiniLAPACK On)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # set the test coverage switch
    set(TestCoverage On)
    set(CMAKE_C_FLAGS_DEBUG "-pg -g3 -O0 -Wall -pedantic -fno-fast-math -fsanitize=address -fsanitize=leak -fopenmp")
    set(CMAKE_Fortran_FLAGS_DEBUG "-pg -g -O0 -Wall -pedantic -Warray-temporaries -fno-omit-frame-pointer -fcheck=all -fbacktrace -fopenmp -finit-real=nan -finit-integer=-9999 -fsanitize=address -fsanitize=leak")
    set(CMAKE_LDFLAGS_DEBUG "-g -fsanitize=address -fsanitize=leak")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
    set(CMAKE_C_FLAGS_DEBUG "-pg -g3 -O0 -Wall -pedantic -fno-fast-math -fsanitize=address -fsanitize=leak -fopenmp")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -Wall -pedantic -fsanitize=address -fno-omit-frame-pointer -fsanitize=leak -fopenmp")
    set(CMAKE_LDFLAGS_DEBUG "-g -fsanitize=address -fsanitize=leak")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
    set(CMAKE_C_FLAGS_DEBUG "-g -O0 -fno-fast-math")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -nan -C=all -C=undefined -u -ieee=full -kind=unique")
    set(CMAKE_LDFLAGS_DEBUG "-g")
  else()
    set(CMAKE_C_FLAGS_DEBUG "-g -O0")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
    set(CMAKE_LDFLAGS_DEBUG "-g")
  endif()
  # set( CMAKE_VERBOSE_MAKEFILE on )
endif()

# Set flags
set(CMAKE_POSITION_INDEPENDENT_CODE ON) # SET fPIC
set(CMAKE_C_STANDARD 99) # SET -std=c99

if(TestCoverage) # Only enabled for GNU compiler
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
   set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
   set(CMAKE_LDFLAGS_DEBUG "${CMAKE_LDFLAGS_DEBUG} -g -fprofile-arcs -ftest-coverage")
   # Provides the new commands: coverage* and clean-coverage
   # INCLUDE(CodeCoverage)
endif(TestCoverage)

# Print build information
message(STATUS "")
message(STATUS "  ${PROJECT_NAME} v${PROJECT_VERSION}")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "")

# Add clean-all target as alias to clean
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
)

# BLAS and LAPACK
if(NOT CompileMiniLAPACK)
  # try to find lapack..
  find_package(LAPACK)
  if(NOT LAPACK_FOUND)
    # couldn't find lapack...compile locally
    set(CompileMiniLAPACK On)
  else()
    message( STATUS "The LAPACK library was found at" ${LAPACK_LIBRARIES} )
    set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
  endif()
endif()

# LAPACK routines
if(CompileMiniLAPACK)
   message("*****************************WARNING**********************************")
   message("* Compiling our own LAPACK library.                                  *")
   message("* This should be avoided as it is likely to be less optimised than a *")
   message("* vendor supplied version, and may also be out of date.              *")
   message("**********************************************************************")
   # note that blas is a dependency of lapack, so blas only explicitly needed here,
   # if we've built our own lapack
   find_package(BLAS REQUIRED)
   if ( NOT BLAS_FOUND )
     message( FATAL_ERROR "BLAS required, but not found!")
   endif ( NOT BLAS_FOUND )
   message( STATUS "The BLAS library was found at" ${BLAS_LIBRARIES} )
   set(LIBS ${LIBS} ${BLAS_LIBRARIES})
   add_subdirectory(lapack)
   set(LAPACK_LIBRARIES ral_nlls_lapack)
   # Set compiler flags for the LAPACK library
   if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
     target_compile_options(ral_nlls_lapack PRIVATE -w -std=legacy)
   elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
     target_compile_options(ral_nlls_lapack PRIVATE -w -std=legacy)
   endif()
   set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
endif()

# Python interface
find_package (Python COMPONENTS Interpreter)

if (PYTHON_Interpreter_FOUND)
   message( STATUS "looking for files in ${CMAKE_CURRENT_SOURCE_DIR}")
   set(PYPROJ_IN "${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml.in")
   set(PYPROJ    "${CMAKE_CURRENT_BINARY_DIR}/pyproject.toml")
   set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
   set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
   set(OUTPUT   "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")

   configure_file(${PYPROJ_IN} ${PYPROJ})
   configure_file(${SETUP_PY_IN} ${SETUP_PY})

   execute_process (
     COMMAND which ${PYTHON_EXECUTABLE}
     OUTPUT_VARIABLE outVar
     )
   message( STATUS ${outVar} )
   message( STATUS "python found -- writing to ${CMAKE_CURRENT_BINARY_DIR}")
   message( STATUS "To install the python library, run: ${PYTHON_EXECUTABLE} -m pip install .")
   add_custom_command(OUTPUT ${OUTPUT}
                     COMMAND ${PYTHON_EXECUTABLE} -m pip install .
                     COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT})
   install(CODE "execute_process(
           COMMAND ${PYTHON_EXECUTABLE} -m pip install .
           )"
          )

   add_custom_target(target ALL DEPENDS ${OUTPUT} ral_nlls)

endif()

# Main libral_nlls library
add_subdirectory(src)

# Tests and examples
add_subdirectory (test)
add_subdirectory (example/C)
add_subdirectory (example/Fortran)
add_subdirectory (example/Python)

if(TestCoverage) # Only enabled for GNU compiler
   INCLUDE(CodeCoverage)
endif(TestCoverage)
# Copyright (c) 2016, The Science and Technology Facilities Council (STFC)
# All rights reserved.
# Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.

# CMake project file for libral_nlls
# libral_nlls Polish: cmake-format --line-width=100
cmake_minimum_required (VERSION 3.25)

project(
  RALFit
  VERSION 2.0.0
  LANGUAGES Fortran C)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

# Option to compile our own LAPACK
option(CompileMiniLAPACK "Compile our own LAPACK (discouraged)" Off)
option(TestCoverage "Perform code coverage analysis" Off)

# Build the library with floating point single precision (WP=float)
option(SINGLE_PRECISION "Build library in SINGLE precision (default Off)" Off)

# Internal option to avoid compiler warning noise issues that wont-be-fixed (unused iface dummies, etc)
set(SUPPRESS_INTERNAL_WONTFIX On)

# Internal option to supress some compiler warning in debug mode (issues from legacy code)
option(SUPPRESS_INTERNAL_FIXMES "Suppress compiler warnings of known issues" Off)

# Address SANitizer option (only for GCC and Clang)
# Debug build will turn it On if option set to "Debug" (default)
set(ASAN "Debug" CACHE STRING "Build with Address / Leak SANitizer (Off/On/Debug)")
set_property(CACHE ASAN PROPERTY STRINGS "Off" "On" "Debug")

# Option to build the Python interface of RALFit
option(BUILD_PYTHON "Build RalNLLS Python wheel (default On)" On)


set(CMAKE_VERBOSE_MAKEFILE On)

if (SINGLE_PRECISION)
    add_compile_definitions(SINGLE_PRECISION)
endif(SINGLE_PRECISION)

# Set the location for cmake module files
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

set(CMAKE_EXE_LINKER_FLAGS "-Wl,-export-dynamic")

set(CMAKE_EXPORT_COMPILE_COMMANDS On)

if(CMAKE_BUILD_TYPE MATCHES Release)
  set(TestCoverage Off)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_C_FLAGS_RELEASE
        "-O2 -Wall -pedantic -march=native -fno-fast-math -fopenmp -fno-omit-frame-pointer")
    set(CMAKE_Fortran_FLAGS_RELEASE
        "-O2 -Wall -pedantic -march=native -fopenmp -fno-omit-frame-pointer")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
    set(CMAKE_C_FLAGS_RELEASE
        "-O2 -Wall -pedantic -march=native -fno-fast-math -fopenmp -fno-omit-frame-pointer")
    set(CMAKE_Fortran_FLAGS_RELEASE
        "-O2 -Wall -pedantic -march=native -fopenmp -fno-omit-frame-pointer")
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
    ##
    set(CMAKE_LDFLAGS_DEBUG "-g3 -fprofile-arcs -ftest-coverage")
    ##
    set(CMAKE_C_FLAGS_DEBUG
        "-pg -g3 -gdwarf-5 -O0 -Wall -pedantic -fno-fast-math -fopenmp")
    set(CMAKE_Fortran_FLAGS_DEBUG
        "-pg -g -O0 -gdwarf-5 -Wall -pedantic -Warray-temporaries -fno-omit-frame-pointer -fcheck=all -fbacktrace -fopenmp -finit-real=nan -finit-integer=-9999"
    )
    set(CMAKE_LDFLAGS_DEBUG "-g")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
    set(CMAKE_C_FLAGS_DEBUG
        "-pg -g3 -O0 -Wall -pedantic -fno-fast-math -fopenmp")
    set(CMAKE_Fortran_FLAGS_DEBUG
        "-g -O0 -Wall -pedantic -fno-omit-frame-pointer -fopenmp"
    )
    set(CMAKE_LDFLAGS_DEBUG "-g")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
    set(CMAKE_C_FLAGS_DEBUG "-g -O0 -fno-fast-math")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -gline -nan -u -ieee=full -kind=unique -C=undefined -C=all")
    set(CMAKE_LDFLAGS_DEBUG "-g")
    # Don't add unsupported debug flags, nagfor compiler does not supports these flags
    set(SUPPRESS_INTERNAL_FIXMES Off)
    set(SUPPRESS_INTERNAL_WONTFIX Off)
  else()
    set(CMAKE_C_FLAGS_DEBUG "-g -O0")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
    set(CMAKE_LDFLAGS_DEBUG "-g")
    # Don't add unsupported debug flags, not known if the compiler supports the flags
    set(SUPPRESS_INTERNAL_FIXMES Off)
    set(SUPPRESS_INTERNAL_WONTFIX Off)
  endif()
endif()

string(TOLOWER "${ASAN}" ASAN)
if((CMAKE_Fortran_COMPILER_ID STREQUAL "Flang") OR (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"))
  if(ASAN STREQUAL "debug" AND CMAKE_BUILD_TYPE MATCHES Debug)
    message(STATUS "Address / Leak SANitizer enabled by Debug build")
    set(ASAN "on")
  endif()
  if(ASAN STREQUAL "on")
    message(STATUS "Address / Leak SANitizer enabled")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fsanitize=leak")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fsanitize=address -fsanitize=leak")
    set(CMAKE_LDFLAGS "${CMAKE_LDFLAGS} -fsanitize=address -fsanitize=leak")
  endif()
endif()

# Set flags
set(CMAKE_POSITION_INDEPENDENT_CODE ON) # SET fPIC
set(CMAKE_C_STANDARD 99) # SET -std=c99
set(CMAKE_Fortran_PREPROCESS ON)

if(TestCoverage) # Only enabled for GNU compiler
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
  set(CMAKE_LDFLAGS_DEBUG "${CMAKE_LDFLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
endif(TestCoverage)

# Print build information
message(STATUS "")
message(STATUS "  ${PROJECT_NAME} v${PROJECT_VERSION}")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
if (SINGLE_PRECISION)
    message(STATUS "  Floating point precision: SINGLE")
else()
    message(STATUS "  Floating point precision: DOUBLE")
endif(SINGLE_PRECISION)
message(STATUS "")

# Add clean-all target as alias to clean
add_custom_target(clean-all COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean)

# BLAS and LAPACK
if(NOT CompileMiniLAPACK)
  # Try to find LAPACK...
  find_package(LAPACK)
  if(NOT LAPACK_FOUND)
    # Couldn't find LAPACK, Compile mini reference BLAS / LAPACK
    set(CompileMiniLAPACK On)
  endif()
endif()

# LAPACK routines
if(CompileMiniLAPACK)
   message("**************************** WARNING *********************************")
   message("* Compiling our mini reference BLAS / LAPACK library.                *")
   message("* This should be avoided as it is likely to be less optimised than a *")
   message("* vendor supplied version.                                           *")
   message("**********************************************************************")
   # Unit-test should only be run/verified when using the provided BLAS and LAPACK
   add_subdirectory(lapack)
   set(LAPACK_LIBRARIES ral_nlls_lapack)
endif()

message(STATUS "The LAPACK / BLAS libraries: " ${LAPACK_LIBRARIES})
set(LIBS ${LIBS} ${LAPACK_LIBRARIES})

# Python interface
if (BUILD_PYTHON)
    find_package(Python COMPONENTS Interpreter)
endif()

if(Python_Interpreter_FOUND)
  message(STATUS "looking for files in ${CMAKE_CURRENT_SOURCE_DIR}")
  set(PYPROJ_IN "${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml.in")
  set(PYPROJ "${CMAKE_CURRENT_BINARY_DIR}/pyproject.toml")
  set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
  set(SETUP_PY "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
  set(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")

  configure_file(${PYPROJ_IN} ${PYPROJ})
  configure_file(${SETUP_PY_IN} ${SETUP_PY})

  message(STATUS "python found -- writing to ${CMAKE_CURRENT_BINARY_DIR}")
  message(STATUS "To install the python library, run: ${Python_EXECUTABLE} -m pip install .")
  add_custom_command(
    OUTPUT ${OUTPUT}
    COMMAND ${Python_EXECUTABLE} -m pip install .
    COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT})
  install(CODE "execute_process(
           COMMAND ${Python_EXECUTABLE} -m pip install .
           )")

  add_custom_target(target ALL DEPENDS ${OUTPUT} ral_nlls)

endif()

# Main libral_nlls library
add_subdirectory(src)

# Tests and examples
add_subdirectory(test)
add_subdirectory(example/C)
add_subdirectory(example/Fortran)
add_subdirectory(example/Python)

if(TestCoverage) # Only enabled for GNU compiler
  include(CodeCoverage)
endif(TestCoverage)

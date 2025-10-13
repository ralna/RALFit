# Copyright (c) 2012 - 2015, Lars Bilke All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted
# provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions
#   and the following disclaimer.
#
# 1. Redistributions in binary form must reproduce the above copyright notice, this list of
#   conditions and the following disclaimer in the documentation and/or other materials provided
#   with the distribution.
#
# 1. Neither the name of the copyright holder nor the names of its contributors may be used to
#   endorse or promote products derived from this software without specific prior written
#   permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# 2012-01-31, Lars Bilke - Enable Code Coverage
#
# 2013-09-17, Joakim Söderberg - Added support for Clang. - Some additional usage instructions.
#
# 2025-10-10, Andrew Sajo - Modernized using new features - Modified to add custom target to run all
# executables for coverage. Polish: cmake-format --line-width=100

# Check prereqs
find_program(LCOV_PATH lcov REQUIRED)
find_program(GENHTML_PATH genhtml REQUIRED)

if(NOT (CMAKE_BUILD_TYPE STREQUAL "Debug"))
  message(FATAL_ERROR "Code coverage on (non-Debug) build is not supported")
endif()

add_custom_target(
  coverage
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMAND ${LCOV_PATH} --directory . --capture --output-file coverage.info
  COMMAND ${LCOV_PATH} --ignore-errors unused --remove coverage.info 'lapack/*' '/usr/*' --exclude '*/lapack/*'
          --output-file coverage.info.cleaned
  COMMAND ${GENHTML_PATH} -o coverage coverage.info.cleaned
  DEPENDS coverage-run-all)

add_custom_target(
  clean-coverage
  COMMAND ${LCOV_PATH} --directory . --zerocounters
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/coverage
  COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/coverage.info
          ${CMAKE_BINARY_DIR}/coverage.info.cleaned
  COMMAND ${CMAKE_COMMAND} -E echo "Resetting coverage counters."
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_dependencies(clean-all clean-coverage)

# Get all executable targets in the project
get_property(
  ALL_TEST_TARGETS
  DIRECTORY "${CMAKE_SOURCE_DIR}/test"
  PROPERTY BUILDSYSTEM_TARGETS)
get_property(
  ALL_EXAMPLE_C_TARGETS
  DIRECTORY "${CMAKE_SOURCE_DIR}/example/C"
  PROPERTY BUILDSYSTEM_TARGETS)
get_property(
  ALL_EXAMPLE_F_TARGETS
  DIRECTORY "${CMAKE_SOURCE_DIR}/example/Fortran"
  PROPERTY BUILDSYSTEM_TARGETS)

# Create custom target to run all executables
add_custom_target(
  coverage-run-all
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMAND ${CMAKE_COMMAND} -E echo "Running all executables for coverage:")

# Add dependencies for all executable targets
foreach(target ${ALL_TEST_TARGETS} ${ALL_EXAMPLE_C_TARGETS} ${ALL_EXAMPLE_F_TARGETS})
  get_target_property(target_type ${target} TYPE)
  if(target_type STREQUAL "EXECUTABLE")
    add_dependencies(coverage-run-all ${target})
    add_custom_command(
      TARGET coverage-run-all
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E echo "Running ${target}..."
      COMMAND $<TARGET_FILE:${target}> > /dev/null 2>&1
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  endif()
endforeach()

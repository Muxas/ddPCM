# @copyright (c) 2020-2020 RWTH Aachen. All rights reserved.
#
# ddX software
#
# @file CMakeLists.txt
# @version 1.0.0
# @author Aleksandr Mikhalev
# @date 2020-12-17


###############################################################################
#                       THIS IS A TOP-LEVEL CMAKELISTS.txt                    #
#                                                                             #
#        It is intended to find all dependencies (required or optional)       #
#                    and set up corresponding variables                       #
###############################################################################


###############################################################################
##                            PRELIMINARIES                                  ##
###############################################################################

# Need to identify lowest possible CMake version
cmake_minimum_required(VERSION 3.2.3)

# Restrict building in top-level directory
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "In-source builds are not allowed.\nPlease create a "
        "build directory first and execute cmake configuration from this "
        "directory. Example: mkdir build && cd build && cmake ..")
endif()

# Read version
file(STRINGS "VERSION.txt" VERSION LIMIT_COUNT 1)

# Notify user about project name and version
message(STATUS "Configuring ddX ${VERSION}")

# Create project and check Fortran compiler
project(ddX VERSION ${VERSION} LANGUAGES Fortran)

###############################################################################
##                        REQUIREMENTS: BLAS AND LAPACK                      ##
###############################################################################

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

###############################################################################
##                           DEFINE OPTIONS                                  ##
###############################################################################

# Option for build type
set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose type of build, options are: None Debug Release\
    RelWithDebInfo MinSizeRel")

# Build tests
option(TESTS "Generate testing binaries" ON)

# Build documentation
option(DOCS "Build documentation" ON)

# Option for a code coverage from cmake/CodeCoverage.cmake
option(CODE_COVERAGE OFF)

###############################################################################
##                               CHECK OPTIONS                               ##
###############################################################################

# Check if testing is required
if(TESTS)
    include_directories(${CMAKE_BINARY_DIR}/src)
    include(CTest)
    enable_testing()
endif()

# Check if need to build documentation
if(DOCS)
    find_package(Doxygen)
    if(NOT DOXYGEN_FOUND)
        set(DOCS OFF)
    endif()
endif()

# Update compiler flags in case of code coverage
if(CODE_COVERAGE)
    set(CMAKE_BUILD_TYPE "Debug")
    # Include code-coverage
    include(${CMAKE_SOURCE_DIR}/cmake/CodeCoverage.cmake)
    append_coverage_compiler_flags()
endif()

###############################################################################
##                 BUILD LIBRARY, EXAMPLES, TESTS AND DOCS                   ##
###############################################################################

add_subdirectory("src")

# Build examples
if(EXAMPLES)
    add_subdirectory("examples")
endif()

# Build tests
if(TESTS)
    add_subdirectory("tests")
endif()

# Build documentation
if(DOCS)
    add_subdirectory("docs")
endif()



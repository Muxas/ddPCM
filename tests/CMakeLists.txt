# @copyright (c) 2020-2020 RWTH Aachen. All rights reserved.
#
# ddX software
#
# @file tests/CMakeLists.txt
# @version 1.0.0
# @author Aleksandr Mikhalev
# @date 2020-02-17

set(TESTS
    "dd_core.f90"
    "main.f90")
list(APPEND COVERAGE_EXCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/*")
set(TESTS_EXE)
message(STATUS "EXCLUDE: ${COVERAGE_EXCLUDES}")

foreach(src ${TESTS})
    # Parse path to get name of executable and its corresponding unique logical
    # name, based on its path and name
    get_filename_component(exec ${src} NAME_WE)
    # Register executable with its unique logical name
    add_executable(${exec} ${src})
    # Link to ddX
    target_link_libraries(${exec} ddx)
    # Register test
    add_test(NAME ${exec} COMMAND ${exec})
    list(APPEND TESTS_EXE ${exec})
endforeach()

# Code coverage
if(CODE_COVERAGE)
    setup_target_for_coverage_lcov(NAME coverage DEPENDENCIES ${TESTS_EXE}
        EXECUTABLE ${exec})
endif()

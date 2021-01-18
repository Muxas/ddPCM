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
    "main.f90"
    "main_ddLPB.f90")
list(APPEND COVERAGE_EXCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/*")
message(STATUS "EXCLUDE: ${COVERAGE_EXCLUDES}")

ddx_add_test("dd_core.f90")
ddx_add_test("main.f90" "${CMAKE_SOURCE_DIR}/tests/Input.txt")
ddx_add_test("main_ddLPB.f90" "${CMAKE_SOURCE_DIR}/tests/Input_ddLPB.txt")


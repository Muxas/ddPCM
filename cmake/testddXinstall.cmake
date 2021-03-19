cmake_minimum_required(VERSION 3.2.3)

project(ddx_test LANGUAGES Fortran)

find_package(ddX REQUIRED)
add_executable(dd_core "dd_core.f90")

set(ddx_incl ${ddX_DIR}/../../include/ddX)
message(${ddx_incl})

include_directories(ddx_test ${ddx_incl})

target_link_libraries(dd_core ddx)
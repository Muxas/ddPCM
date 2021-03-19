cmake_minimum_required(VERSION 3.2.3)

project(ddx_test LANGUAGES Fortran)

find_package(ddX REQUIRED)
add_executable(dd_core "dd_core.f90")

target_link_libraries(dd_core ddx)
# @copyright (c) 2020-2020 RWTH Aachen. All rights reserved.
#
# ddX software
#
# @file src/CMakeLists.txt
# @version 1.0.0
# @author Aleksandr Mikhalev
# @date 2020-02-17

set(SRC
    "dd_core.f90"
    "dd_operators.f90"
    "llgnew.f")

set(COVERAGE_EXCLUDES
    "${CMAKE_CURRENT_SOURCE_DIR}/llgnew.f" PARENT_SCOPE)

list(TRANSFORM SRC PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/)

# Define ddX library
add_library(ddx SHARED ${SRC})
# Link against BLAS and LAPACK
target_link_libraries(ddx PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
# Code coverage
if(CODE_COVERAGE)
    target_link_libraries(ddx PRIVATE gcov)
endif()

# Move SRC to parent directory
set(SRC ${SRC} PARENT_SCOPE)


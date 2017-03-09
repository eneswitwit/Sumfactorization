cmake_minimum_required (VERSION 2.6)
project (Sumfactorization)

# set the output files location
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

# add the binary tree to the search path for include files
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SORUCE_DIR}/include")

# add executable
add_executable(testvmult ${PROJECT_SOURCE_DIR}/tests/integrate.cc)





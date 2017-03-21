cmake_minimum_required (VERSION 2.6)
project (Sumfactorization)

# set the output files location
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

# add the binary tree to the search path for include files
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SORUCE_DIR}/include")

# check for a sufficiently recent compiler
INCLUDE(CheckCXXCompilerFlag)

CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
CHECK_CXX_COMPILER_FLAG("-std=c++1y" COMPILER_SUPPORTS_CXX1Y)
IF(COMPILER_SUPPORTS_CXX14)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
ELSEIF(COMPILER_SUPPORTS_CXX1Y)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")
ELSE()
  MESSAGE(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
ENDIF()

# use reasonable flags for warnings
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# add executable
add_executable(testvmult ${PROJECT_SOURCE_DIR}/tests/integrate.cc)
add_executable(polynomial ${PROJECT_SOURCE_DIR}/tests/polynomial.cc)
add_executable(quadrature ${PROJECT_SOURCE_DIR}/tests/quadrature.cc)

# Put all tests in the test directory, where the sources also are
set(RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)

# testing
enable_testing()

# check if test runs and exits without error
add_test(testvmult ${PROJECT_SOURCE_DIR}/bin/testvmult)
add_test(polynomial ${PROJECT_SOURCE_DIR}/bin/polynomial)
add_test(quadrature ${PROJECT_SOURCE_DIR}/bin/quadrature)





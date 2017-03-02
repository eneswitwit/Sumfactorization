cmake_minimum_required (VERSION 2.6)
project (integrate)

#FILE(GLOB sources source/*.cc)
INCLUDE_DIRECTORIES(.)

add_executable(integrate Integrate.cpp)
#add_executable(test_explicit test_explicit.cc ${sources})



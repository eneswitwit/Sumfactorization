cmake_minimum_required (VERSION 2.6)
project (Integrate)

# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")

# enable dashboard scripting
include (CTest)

# add the executable
add_executable(Integrate integrate.cpp)

# enable testing
enable_testing ()

# does the application run
add_test (IntegrationRuns Integrate)

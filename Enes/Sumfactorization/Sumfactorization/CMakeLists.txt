cmake_minimum_required (VERSION 2.6)
project (Sumfactorization)

# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")

# add the executable
add_executable(Sumfactorization tests/integrate.cc)

# target Libraries
# target_link_libraries (Sumfactorization include)

# Install rules
# install (TARGETS Sumfactorization DESTINATION bin)

# enable dashboard scripting
include (CTest)

# enable testing
enable_testing ()

# does the application run
add_test (SumfactorizationRuns Sumfactorization)
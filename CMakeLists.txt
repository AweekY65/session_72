cmake_minimum_required(VERSION 3.16)
project(local_octree CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-Wall -Wextra -O2)

add_library(octree src/octree.cpp)
target_include_directories(octree PUBLIC include)

add_executable(test_octree tests/test_octree.cpp)
target_link_libraries(test_octree PRIVATE octree)

enable_testing()
add_test(NAME octree_tests COMMAND test_octree)

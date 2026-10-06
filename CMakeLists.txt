cmake_minimum_required(VERSION 3.16)
project(local_octree CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_library(octree src/octree.cpp)
target_include_directories(octree PUBLIC include)
target_compile_options(octree PRIVATE -Wall -Wextra)

enable_testing()
add_executable(test_octree tests/test_octree.cpp)
target_link_libraries(test_octree PRIVATE octree)
add_test(NAME octree_tests COMMAND test_octree)

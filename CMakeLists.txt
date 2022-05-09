cmake_minimum_required(VERSION 3.22)
project(untitled)

set(CMAKE_CXX_STANDARD 20)
add_subdirectory(lib)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(untitled main.cpp Calculation.cpp Calculation.h)
target_link_libraries(untitled gtest gtest_main)
cmake_minimum_required(VERSION 3.15)
project(untitled2)

SET(CMAKE_CXX_COMPILER ../../afl-clang-fast++)

add_executable(simple_test main.cpp)
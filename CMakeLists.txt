cmake_minimum_required(VERSION 3.20)
project(Lab_7_vectorize C)

set(CMAKE_C_STANDARD 99)

set(CMAKE_C_FLAGS "-g -O3 -msse3")
add_executable(Lab_7_vectorize  without.c)

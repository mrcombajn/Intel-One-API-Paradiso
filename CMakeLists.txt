cmake_minimum_required(VERSION 4.2)
project(paradiso)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(paradiso
        minimal.cpp
        minimal.h)

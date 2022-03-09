cmake_minimum_required(VERSION 3.5)
project(10)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
        life2d.c)

add_executable(10 ${SOURCE_FILES})
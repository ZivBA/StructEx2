cmake_minimum_required(VERSION 3.6)
project(StructEx2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(./gamb)
file(GLOB_RECURSE SOURCE_FILES main.cpp "gamb/*.h" "gamb/*.cc")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "/home/rooty/ClionProjects/StructEx2/input_files")
#set(SOURCE_FILES main.cpp )
add_executable(StructEx2 ${SOURCE_FILES})
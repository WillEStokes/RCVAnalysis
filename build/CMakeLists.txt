# cmake system name
set(CMAKE_SYSTEM_PROCESSOR x86_64)

# set paths as variable
set(toolchainpath C:/MinGW)

# configure compilers
set(CMAKE_CXX_COMPILER ${toolchainpath}/bin/g++.exe)

# compiler flags
set(CMAKE_<LANG>_STANDARD_INCLUDE_DIRECTORIES ${toolchainpath}/lib/gcc/mingw32/9.2.0/include ${toolchainpath}/lib/gcc/mingw32/9.2.0/include/c++ ${toolchainpath}/bin)

# specify search behaviour
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

# linker options
set(CMAKE_CXX_STANDARD_LIBRARIES "-static-libgcc -static-libstdc++")

# project specific information
cmake_minimum_required(VERSION 3.17)
project(arrayPractice)
# set(CMAKE_CXX_FLAGS "-std=c++0x -Wall")
# set(CMAKE_EXE_LINKER_FLAGS "-static")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_BUILD_TYPE Debug)
# add_library(arrayPractice SHARED
add_executable(arrayPractice WIN32
../src/arrayPractice.cpp)
# ../src/arrayPractice.h)
# target_sources(arrayPractice PRIVATE)
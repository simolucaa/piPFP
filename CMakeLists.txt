cmake_minimum_required(VERSION 3.10)

# Set the project name
project(piPFP)

add_subdirectory(growt)

# Add the executable
add_executable(piPFP piPFP.cpp piPFP.hpp utils.c xerrors.c)

add_executable(piPFP_growth piPFP_growth.cpp piPFP.hpp utils.c xerrors.c)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -march=native -g")

# include the growt library
target_include_directories(piPFP PUBLIC growt)
target_include_directories(piPFP_growth PUBLIC growt)

find_package(Threads REQUIRED)
target_link_libraries(piPFP pthread)
target_link_libraries(piPFP_growth pthread)


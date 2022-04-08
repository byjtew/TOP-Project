cmake_minimum_required(VERSION 3.16)

# Project information
project(top-project LANGUAGES C)
# Todo: add project description
set(PROJECT_DESCRIPTION "")

## Build type (Debug, Release, RelWithDebInfo, MinSizeRel)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

## C standards
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

## C flags
set(CMAKE_C_FLAGS "-Wall -Wextra -pedantic -std=c99 -lm")

## MPI support
message(STATUS "Looking for mpicc...")
if (NOT CMAKE_C_COMPILER MATCHES ".*mpicc")
    message(FATAL_ERROR "CMAKE_C_COMPILER must be 'mpicc'\nRebuild using: -DCMAKE_C_COMPILER=`which mpicc`")
else ()
    message(STATUS "Found mpicc: ${CMAKE_C_COMPILER}")
endif ()

## OpenMP support
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message(STATUS "OpenMP found")
else ()
    message(STATUS "OpenMP not found")
endif ()
find_package(Threads REQUIRED)
if (Threads_FOUND)
    message(STATUS "Threads found")
else ()
    message(STATUS "Threads not found")
endif ()

# Libraries
message(STATUS "Looking for libraries...")
set(INC_DIR ${CMAKE_SOURCE_DIR}/include)
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
add_subdirectory(src)

# Executables
set(CONFIG_FILE ${CMAKE_SOURCE_DIR}/config.txt)
add_subdirectory(bin)

# Testing
option(TEST "Enable testing" OFF)
if (TEST)
    message(STATUS "Enabling testing")
    enable_testing()
    add_subdirectory(tests)
else ()
    message(STATUS "Testing disabled. To enable, use: -DTEST=ON")
endif ()
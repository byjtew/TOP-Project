cmake_minimum_required(VERSION 3.16)

# Project information
project(top-project LANGUAGES C)
# Todo: add project description
set(PROJECT_DESCRIPTION "")

## Build type (Debug, Release, RelWithDebInfo, MinSizeRel)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif ()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

## C standards
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

## MPI support
message(STATUS "Looking for mpicc...")
if (NOT CMAKE_C_COMPILER MATCHES ".*mpicc")
    message(FATAL_ERROR "CMAKE_C_COMPILER must be 'mpicc'\nRebuild using: -DCMAKE_C_COMPILER=`which mpicc`")
else ()
    message(STATUS "Found mpicc: ${CMAKE_C_COMPILER}")
endif ()

## OpenMP support
find_package(Threads REQUIRED)
if (Threads_FOUND)
    message(STATUS "Threads found")
else ()
    message(FATAL_ERROR "Threads not found")
endif ()
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message(STATUS "OpenMP found")
else ()
    message(MESSAGE "OpenMP not found")
endif ()
if (OPENMP_C_FOUND)
    message(STATUS "OpenMP_C found")
else ()
    message(FATAL_ERROR "OpenMP_C not found")
endif ()

## C flags
if (${CMAKE_BUILD_TYPE} MATCHES Debug)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -g -O0 -Wall -Wextra -Wpedantic")
else ()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O3 -ffast-math -march=native -Wall -Wextra -Wpedantic")
endif ()
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
message(STATUS "C flags: ${CMAKE_C_FLAGS}")

# Libraries
message(STATUS "Looking for libraries...")
set(INC_DIR ${CMAKE_SOURCE_DIR}/include)
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
add_subdirectory(src)

# Executables
configure_file(${CMAKE_SOURCE_DIR}/config.txt ${CMAKE_BINARY_DIR}/bin/config.txt COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/small_config.txt ${CMAKE_BINARY_DIR}/bin/small_config.txt COPYONLY)
add_subdirectory(bin)

# Scripts
configure_file(${CMAKE_SOURCE_DIR}/gen_animate_gif_legacy.sh ${CMAKE_BINARY_DIR}/bin/gen_animate_gif_legacy.sh COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/gen_animate_gif.sh ${CMAKE_BINARY_DIR}/bin/gen_animate_gif.sh COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/gen_single_image.sh ${CMAKE_BINARY_DIR}/bin/gen_single_image.sh COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/parallel_gen_animate_gif_legacy.sh ${CMAKE_BINARY_DIR}/bin/parallel_gen_animate_gif_legacy.sh COPYONLY)

# Testing
option(TEST "Enable testing" OFF)
if (TEST)
    message(STATUS "Enabling testing")
    enable_testing()
    add_subdirectory(tests)
else ()
    message(STATUS "Testing disabled. To enable, use: -DTEST=ON")
endif ()
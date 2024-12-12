# CMakeLists.txt based on the Makefile of the old version and CMakeLists.txt of cavity2d
# of palabos.

cmake_minimum_required(VERSION 3.10)

project(embolism)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS True)

include(helper.cmake)

# Set the path for dependency libraries
#
# Note: The name LAMMPS_ROOT will be ignored by CMake due to "policy CMP0074". This is a
# policy regarding how CMake handles package-specific root variables. So I prefixed all
# the names with "PROJ_"
set_env_var_or_default(PROJ_PALABOS_ROOT
    "$ENV{HOME}/git/projects/palabos/bugfix-pointer-access-offLatticeBoundaryCondition3D")
set_env_var_or_default(PROJ_LAMMPS_ROOT
    "$ENV{HOME}/git/projects/lammps/lammps-ext-fix_fcm")
set_env_var_or_default(PROJ_PALABOS_LAMMPS_ROOT
    "$ENV{HOME}/git/projects/palabos-lammps/dev")

# This seems to be OK because unlike LAMMPS, apparently we are not looking for Palabos
# package in the current source file.
set(PALABOS_ROOT ${PROJ_PALABOS_ROOT})

# # Path to palabos
# # set(PROJ_PALABOS_ROOT "$ENV{HOME}/git/software/lbm-md/palabos")
# if(DEFINED ENV{PROJ_PALABOS_ROOT})
#     set(PROJ_PALABOS_ROOT "$ENV{PROJ_PALABOS_ROOT}")
# else()
#     set(PROJ_PALABOS_ROOT "$ENV{HOME}/git/projects/palabos/bugfix-pointer-access-offLatticeBoundaryCondition3D")
# endif()
# message(STATUS "PROJ_PALABOS_ROOT is set to ${PROJ_PALABOS_ROOT}")
# set(PALABOS_ROOT ${PROJ_PALABOS_ROOT})
#
# # Path to lammps
# # set(PROJ_LAMMPS_ROOT "$ENV{HOME}/git/software/lbm-md/lammps")
# if(DEFINED ENV{PROJ_LAMMPS_ROOT})
#     set(PROJ_LAMMPS_ROOT "$ENV{PROJ_LAMMPS_ROOT}")
# else()
#     set(PROJ_LAMMPS_ROOT "$ENV{HOME}/git/projects/lammps/lammps-ext-fix_fcm")
# endif()
# message(STATUS "PROJ_LAMMPS_ROOT is set to ${PROJ_LAMMPS_ROOT}")
#
# # Path to palabos-lammps
# # set(PROJ_LAMMPS_ROOT "$ENV{HOME}/git/software/lbm-md/palabos-lammps")
# if(DEFINED ENV{PROJ_PALABOS_LAMMPS_ROOT})
#     set(PROJ_PALABOS_LAMMPS_ROOT "$ENV{PROJ_PALABOS_LAMMPS_ROOT}")
# else()
#     set(PROJ_PALABOS_LAMMPS_ROOT "$ENV{HOME}/git/projects/palabos-lammps/dev")
# endif()
# message(STATUS "PROJ_PALABOS_LAMMPS_ROOT is set to ${PROJ_PALABOS_LAMMPS_ROOT}")

# Include directories
include_directories(
    ${PROJ_PALABOS_ROOT}/src
    ${PROJ_PALABOS_ROOT}/externalLibraries
    ${PROJ_LAMMPS_ROOT}/src
    ${PROJ_PALABOS_LAMMPS_ROOT}/src
)

# Link libraries
link_directories(
    $ENV{HOME}/.local/lib
)

# Add the LAMMPS CMake directory to the search path
set(CMAKE_PREFIX_PATH "/Users/ali/.local/lib/cmake" ${CMAKE_PREFIX_PATH})

# Find the LAMMPS package
find_package(LAMMPS REQUIRED)

# Source files for this project
set(PROJECT_FILES embolism.cpp)

# Compiler options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wnon-virtual-dtor")
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

# MPI support
option(MPI_PARALLEL "Enable MPI parallelism" ON)
if(MPI_PARALLEL)
    find_package(MPI REQUIRED)
    include_directories(${MPI_INCLUDE_PATH})
    set(CMAKE_CXX_COMPILER mpicxx)
    add_definitions(-DPLB_MPI_PARALLEL)
endif()

# POSIX support
option(ENABLE_POSIX "Enable POSIX API" ON)
if(ENABLE_POSIX)
    add_definitions(-DPLB_USE_POSIX)
endif()

# Apple-specific definitions
if(APPLE)
    add_definitions(-DPLB_MAC_OS_X)
endif()

# Windows-specific definitions
if(WIN32)
    add_definitions(-DPLB_WINDOWS)
endif()

# Palabos library
file(GLOB_RECURSE PALABOS_SRC "${PALABOS_ROOT}/src/*.cpp")
file(GLOB_RECURSE PALABOS_EXTERNAL_SRC "${PALABOS_ROOT}/externalLibraries/tinyxml/*.cpp")
add_library(palabos STATIC ${PALABOS_SRC} ${PALABOS_EXTERNAL_SRC})

# Project executable
add_executable(embolism ${PROJECT_FILES})

# Link with Palabos and external libraries
target_link_libraries(embolism palabos LAMMPS::lammps)
if(MPI_PARALLEL)
    target_link_libraries(embolism ${MPI_LIBRARIES})
endif()


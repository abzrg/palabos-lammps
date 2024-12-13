cmake_minimum_required(VERSION 3.10)

project(embolism)

# ---------------------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---------------------------------------------------------------------------------------

set(PALABOS_LAMMPS_ROOT "../..")

# ---------------------------------------------------------------------------------------
# Dependencies

# Plabos target
include(${PALABOS_LAMMPS_ROOT}/cmake/palabos.cmake)

# Lammps interface target
include(${PALABOS_LAMMPS_ROOT}/cmake/lammps.cmake)

# Lammps extension target
add_subdirectory(
    ${PALABOS_LAMMPS_ROOT}/lammps-ext ${CMAKE_CURRENT_BINARY_DIR}/lammps-ext
)

# ---------------------------------------------------------------------------------------

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

# OS-specific definitions
if(APPLE)
    add_definitions(-DPLB_MAC_OS_X)
elseif(WIN32)
    add_definitions(-DPLB_WINDOWS)
endif()

# ---------------------------------------------------------------------------------------

add_executable(embolism embolism.cpp)

# Coupling sources
include_directories(${PALABOS_LAMMPS_ROOT}/src)

target_link_libraries(embolism
    PRIVATE palabos
            LAMMPS::lammps
            lammps_ext)

if(MPI_PARALLEL)
    target_link_libraries(embolism
        PRIVATE ${MPI_LIBRARIES})
endif()


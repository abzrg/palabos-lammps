cmake_minimum_required(VERSION 3.10)

project(lammps_ext)

include(../cmake/lammps.cmake)  # target: lammps

add_library(lammps_ext STATIC)

target_sources(lammps_ext
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fix_fcm.cpp)
target_include_directories(lammps_ext
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(lammps_ext PUBLIC LAMMPS::lammps)


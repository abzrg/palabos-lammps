set(HOME $ENV{HOME})

include(${CMAKE_CURRENT_LIST_DIR}/lammps.cmake)

if(NOT TARGET palabos)
    add_library(palabos STATIC)

    # List of required sources files
    file(GLOB_RECURSE PALABOS_SRC "$ENV{PALABOS_ROOT}/src/*.cpp")
    file(GLOB_RECURSE EXT_SRC "$ENV{PALABOS_ROOT}/externalLibraries/tinyxml/*.cpp")

    target_sources(palabos
        PRIVATE ${PALABOS_SRC} ${EXT_SRC})

    target_include_directories(palabos
        PUBLIC $ENV{PALABOS_ROOT}/src
               $ENV{PALABOS_ROOT}/externalLibraries)

    target_link_libraries(palabos
        PUBLIC LAMMPS::lammps)
endif()

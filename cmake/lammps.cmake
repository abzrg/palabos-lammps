set(HOME $ENV{HOME})

if(NOT TARGET lammps_interface)
    add_library(lammps_interface INTERFACE)
    # The following is based on the output of the command:
    # `pkg-config --cflags --libs liblammps`
    target_include_directories(lammps_interface INTERFACE "$ENV{LAMMPS_INSTALL_PREFIX}/include")
    target_compile_definitions(lammps_interface INTERFACE "LAMMPS_SMALLBIG")
    target_link_directories(lammps_interface INTERFACE    "$ENV{LAMMPS_INSTALL_PREFIX}/lib")
    target_link_libraries(lammps_interface INTERFACE lammps)

    # lammps_interface is tempoary and will be replaced by targets provided by the LAMMPS
    # itself through findpackage. For unknown reasons the current targets are not
    # functioning correctly.
    add_library(LAMMPS::lammps ALIAS lammps_interface)
endif()

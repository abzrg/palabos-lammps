# ---------------------------------------------------------------------------------------

# EDIT HERE
# - Please edit the following variables and provide the correct required path.
# - Alternatively set corresponding environment variables in the shell.

LAMMPS_INSTALL_PREFIX := $$HOME/.local
PALABOS_ROOT := \
	$$HOME/git/projects/palabos/bugfix-pointer-access-offLatticeBoundaryCondition3D

# ---------------------------------------------------------------------------------------

ENVS := LAMMPS_INSTALL_PREFIX=$(LAMMPS_INSTALL_PREFIX) PALABOS_ROOT=$(PALABOS_ROOT)
CMAKE := $(ENVS) cmake
BUILDDIR := build
__CMAKE_CACHE := $(BUILDDIR)/CMakeCache.txt

$(if $(filter __%, $(MAKECMDGOALS)), \
	$(error targets prefixed with '__' are only for internal use))

# ---------------------------------------------------------------------------------------

all:
	$(CMAKE) --build $(BUILDDIR)


$(_CMAKE_CACHE): CMakeLists.txt
	$(CMAKE) -B $(BUILDDIR) -S .


clean:
	cmake --build $(BUILDDIR) --target $@


.PHONY: all clean

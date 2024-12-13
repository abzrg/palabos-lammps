# Embolism

## Getting Started

### Environment Variables

`embolism` depends on Palabos, LAMMPS, and a LAMMPS extension and a few coupling code in this repository. For the first two one must specify the two paths in the `Makefile`.
These are:

1. `LAMMPS_INSTALL_PREFIX`: Path to the installation of the LAMMPS library
2. `PALABOS_ROOT`: Path to the root of the Palabos source code.

For example:

```make
# file: Makefile

LAMMPS_INSTALL_PREFIX := $$HOME/.local
PALABOS_ROOT := $$HOME/projects/palabos
```

Alternatively, one can set corresponding environment variables in the shell.

```sh
# In command line

export LAMMPS_INSTALL_PREFIX="$HOME/.local"
export PALABOS_ROOT="$HOME/projects/palabos"
```


### Build and Run

```sh
# Build
make --parallel <N>  # <N>: number of parallel jobs

# Run
mpirun -np <N> build/embolism in.embolism
```


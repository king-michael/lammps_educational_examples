#!/bin/bash

# # serial
# lmp -in input.lammps

# # with mpirun
# mpirun lmp -in input.lammps

# with mpirun and OMP styles (faster)
mpirun lmp -in input.lammps -sf omp -pk omp 1
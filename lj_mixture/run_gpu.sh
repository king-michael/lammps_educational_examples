#!/bin/bash

# mpirun with gpu
# mpirun lmp -in input.lammps -sf gpu -pk gpu 1

lmp -in input.lammps -sf gpu -pk gpu 1
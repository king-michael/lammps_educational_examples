#!/bin/bash

# # serial
#        python script  LAMMPS script couple_frequency
# python vizplotgui_pymol.py input.lammps 100

# parallel
#              python script  LAMMPS script couple_frequency
mpirun python vizplotgui_pymol.py input.lammps 100

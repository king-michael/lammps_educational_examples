# Details

LJ A - B mixture.

liquid at temperature =  0.741


# HowTo
Examples how to run the script

## MPI

$ ./run.sh  # predefine command
$ mpirun lmp -in input.lammps # mpirun
$ mpirun lmp -in input.lammps -sf omp -pk omp 1 # mpirun + accelerated OMP styles

## GPU

$ ./run_gpu.sh # predefine command
$ lmp -in input.lammps -sf gpu -pk gpu 1 # single core + gpu (faster then MPI!)


#Interaction Counter

## Purpose ##
This example creates a set of particles, and counts how many neighbours each particle has.
This can be used for verification/testing of different implementations in the future, 
and is initially being used to show the way a user might create a program using the DSL.

## Defaults ##
This sets up the environment required for the example. To use it, copy it to the root
directory of the repository.

## Initialisation ##
At the moment, the initialisation is done in the `interaction_count_init.rg` file. 
This initialises 9 particles in a grid, with a massive cutoff (such that all are 
in range of each other).
This will be updated to use the HDF5 file reader once implemented, and create a 
variety of tests to check neighbour search is functioning correctly.

## Particle type ##
The particle type is defined in the `interaction_count_part.rg` file. It uses a
standard particle type with one additional field, the interactions field, used to 
store how many interactions each particle has.

## Kernels ##
The kernels are defined in the `interaction_count_kernel.rg` file. There is a 
symmetric kernel, and an asymmetric kernel that is not currently used.

## Main program ##
The main program is specified in the `interaction_count_program.rg` file.
Before the `main_task`, we create the tasks used for the program using the appropriate
functions. In the `main_task`, we create the particles (which the IO system or other
yet to be created functionality will do instead), do a small amount of housekeeping, and
then launch the tasks. We use an assert to check correctness for now.

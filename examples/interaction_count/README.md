# Interaction Counter

## Purpose ##
This example creates a set of particles, and counts how many neighbours each particle has.
This can be used for verification/testing of different implementations in the future, 
and is initially being used to show the way a user might create a program using the DSL.

## Defaults ##
This sets up the environment required for the example. To use it, copy it to the root
directory of the repository.

# Infrastructure

## Info ##
The `infrastructure` directory contains code that will eventually be part of the DSL, such as 
hdf5 file management, particle initialisation etc.

## Initialisation ##
At the moment, the initialisation is done in the `interaction_count_init.rg` file. 
This initialises 9 particles in a grid, with a large cutoff (such that all are 
in range of each other).

## Variables ##
The `interaction_count_variables.rg` file contains some variable declarations, which will be used at a later
date to enable metaprogramming to remove most of the boilerplate code in the main program

# User code
The main directory contains the code that will be required of the user for custom particle types or custom
kernels, as well as the main program used.

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
functions. In the `main_task`, we create the particles using an IO module, and
then launch the tasks. We use an assert to check correctness for now, and then do some HDF5 testing. The correctness checking and HDF5 business
will be moved into the infrastructure at a later date.

The only rule for the main task itself (or anywhere else the user code requires `[function]` style syntax) is that each line must end in a `;` 
to ensure Regent can correctly parse the code.

The main program uses the `simple_HDF5` IO module to initialise the data structures.

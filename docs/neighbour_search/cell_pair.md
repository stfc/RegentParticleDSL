# Cell Pair Algorithm

The cell pair algorithm is an implementation of the cell-lists algorithm for CPU computation. The algorithm is 
periodic, and is found in the `src/neighbour_search/*_cell_pair_tradequeues_*` directories.

## Required Headers
The headers are automatically imported through the `setup_part` and `setup_dsl` function calls used to initialise the framework, according
to the periodicity and dimensionality requirements requested. 

## Algorithm initialisation and maintenance
To initalise the neighbour search algorithm, the following code is required in the main of the program (after the IO initialisation is completed):
```
  [neighbour_init.initialise(variables)];
  [neighbour_init.update_cells(variables)];
```

## The Invoke Syntax

Instead of manually creating and running tasks, HartreeParticleDSL now uses a new type of  function, called Invoke. To use the Invoke
functionality, you pass the config, some tuples containing the kernel function and a type descriptor, followed by an optional
argument to ask for a barrier at the end of the invoke call. For example, if we have two symmetric pairwise functions (named
`sym1` and `sym2`) and a per-particle timestepping kernel (named `timestep`) we would do:
```
[invoke(variable.config, {sym1, SYMMETRIC_PAIRWISE}, {sym2, SYMMETRIC_PAIRWISE}, {timestep, PER_PART}, NO_BARRIER)];
```

When reached, this code will analyse the kernels, and combine kernels into the same task where it deems possible. The criteria for 
combining kernels into the same task are:
1. No Read after Write, Write after Write or Read after Write dependency.
2. The kernel will not result in the neighbour search being invalidated (i.e. the particles are not moved, and the cutoff radius is
   not modified).

For per-particle tasks, these criteria can be ignored, but at the moment they are treated no differently to other task types.

These tasks are combined in the order the kernels appear in the invoke. Any kernel that doesn't meet either criteria for combining
kernels are run in a stand-alone task, though in theory kernels that only break the first rule can still be combined with following 
kernels.

To maximize the combining of kernels, any kernels that do not require a specific ordering should be grouped by type in the invoke call.

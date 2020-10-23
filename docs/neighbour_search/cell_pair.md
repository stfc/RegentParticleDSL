# Cell Pair Algorithm

The cell pair algorithm is an implementation of the cell-lists algorithm for CPU computation. The algorithm is 
periodic, and is found in the `src/neighbour_search/cell_pair_tradequeues` directory.

## Required Headers
The required import to use the algorithm is
```
require("src/neighbour_search/cell_pair_tradequeues/import_cell_pair")
local neighbour_init = require("src/neighbour_search/cell_pair_tradequeues/neighbour_init")
```

## Algorithm initialisation and maintenance
The neighbour search algorithms don't yet have an initialisation or variable. To initalise the neighbour search algorithm,
the following code is required in the main of the program (after the IO initialisation is completed):
```
  [neighbour_init.initialise(variables)];
  [neighbour_init.update_cells(variables)];
```

## Kernel options
This implementation provides support for 3 types of kernel at current.

### Per-particle Kernel
The per-particle Kernel supported by the cell pair algorithm has the following signature:
```
function kick_kernel(part, config)
  local kernel = rquote
    --Kernel code goes here
  end
  return kernel
end
```
The `part` is the `part` fspace type declared by the user, and the `config` argument is the `config` fspace type.

To declare the launch of a set of per-particle tasks, the algorithm provides the following function:
```
local per_part_task = run_per_particle_task( kernel_name, variables.config, neighbour_init.cell_partition )
```

To run the created tasks, use the following signature:
```
  [per_part_task];
```

### Symmetric pairwise task
The symmetric pairwise task has the following signature:
```
function symmetric_pairwise_kernel(part1, part2, r2)
  local kernel = rquote
    --Kernel code goes here
  end
  return kernel
end
```

Both `part1` and `part2` are both `part` fspace types, and both are fully modifiable. 
`r2` is a floating point value, and contains the square of the distance between the particles. Changes to r2 are not visible outside the kernel.

The cutoff radius used for symmetric pairwise tasks is the maximum of the cutoff radii of `part1` and `part2`.

To create the code to run a set of symmetric tasks, the algorithm provides the following function:
```
local symmetric_task = create_symmetric_pairwise_runner( kernel_name, variables.config, neighbour_init.cell_partition )
```

To run the created task, use the following signature:
```
  [symmetric_task];
```

### Asymmetric pairwise task
The asymmetric pairwise task has the following signature:
```
function asymmetric_pairwise_kernel(part1, part2, r2)
  local kernel = rquote
    --Kernel code goes here
  end
  return kernel
end
```

Both `part1` and `part2` are both `part` fspace types. Modifications are only allowed to `part1`
`r2` is a floating point value, and contains the square of the distance between the particles. Changes to r2 are not visible outside the kernel.

The cutoff radius used for asymmetric pairwise tasks is the cutoff radius of `part1`.

To create a asymmetric task, the algorithm provides the following function:
```
local asymmetric_task = create_asymmetric_pairwise_runner( kernel_name, variables.config, neighbour_init.cell_partition )
```

To run the created task, use the following signature:
```
  [asymmetric_task];
```

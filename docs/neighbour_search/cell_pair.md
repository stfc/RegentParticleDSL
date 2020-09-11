# Cell Pair Algorithm

The cell pair algorithm is an implementation of the cell-lists algorithm for CPU computation. The algorithm is 
periodic only, and is found in the `src/neighbour_search/cell_pair` directory.

## Required Headers
The required import to use the algorithm is
```
require("src/neighbour_search/cell_pair/import_cell_pair")
```

## Algorithm initialisation and maintenance
The neighbour search algorithms don't yet have an initialisation or variable. To initalise the neighbour search algorithm,
the following code is required in the main of the program (after the IO initialisation is completed):
```
initialise_cells(variables.config , variables.particle_array)
particles_to_cell_launcher( variables.particle_array, variables.config)

var cell_partition1 = update_cell_partitions([variables.particle_array], [variables.config])
```

Upon completion of any task that moves the particles, the cell lists need to be update with this call:
```
__delete(cell_partition1)
var cell_partition2 = update_cell_partitions(variables.particle_array, variables.config)
```

Note that the same variable cannot be used in the same scope for different `update_cell_partitions` calls. Instead use `__delete` to
ensure safe cleanup, and create a new variable to store the result of the call.

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

To create a per-particle task, the algorithm provides the following function:
```
local per_part_task = generate_per_part_task( kernel_name )
```

To run the created task, use the following signature:
```
  per_part_task(variables.particle_array, cell_partition, variables.config)
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

To create a symmetric task, the algorithm provides the following function:
```
local symmetric_task = create_symmetric_pairwise_runner( kernel_name )
```

To run the created task, use the following signature:
```
  symmetric_task(variables.particle_array, cell_partition, variables.config)
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
local asymmetric_task = create_asymmetric_pairwise_runner( kernel_name )
```

To run the created task, use the following signature:
```
  asymmetric_task(variables.particle_array, cell_partition, variables.config)
```

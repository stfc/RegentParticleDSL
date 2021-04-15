# DSL Program Structure

HartreeParticleDSL programs consist of three sections:

1. The particle type declaration
2. The kernel declaration(s)
3. The main program declaration.

Additional functionality can be added to programs, but these are the minimum required.

## Importing the DSL and setup.
When using HartreeParticleDSL, the system is imported using the following header in the main program file:
```
  import "regent"
  require("src/RegentParticleDSL")
```
Once the header is imported, the system provides a set of calls for setting up the particle method. These are:
- `set_dimensionality(int)`: This call is used to set the dimensionality of the system. Currently the supported values are
`2` or `3`. By default this is set to `3`.
- `set_periodicity(bool)`: This call is used to enable or disable periodicity. The default is `true`.
- `enable_timing(bool)`: This call is used to enable or disable the timing infrastructure. The default is `false`. Note: 
Enabling the timing infrastructure can currently cause additional synchronicity, so should only be used if significant
performance issues are occurring.

Once these functions have been called with appropriate values, the user must call the `setup_part()` function. This prepares the
DSL's internal structures to enable declaration of the particle type and kernel functions.

Once the particle type has been declared, the final setup call required is to `setup_dsl()`. This imports the remaining headers required
by the DSL, and allows the main function to be created.

The overall setup for the DSL should therefore be something like:
```
  import "regent"
  require("src/RegentParticleDSL")
  set_dimensionality(2)
  set_periodicity(true)
  enable_timing(false)
  setup_part()

  --Declare/import particle type
  setup_dsl()

...
``` 

## The Particle Type

Regent's particle type is the main data structure used in the particle methods. A base declaration is
available in `src/particles/default_part.rg`. The outline of the particle type is:
```
  import "regent"
  
  fspace part{
    neighbour_part_space : neighbour_part,
    core_part_space : core_part
  }
```

The `fspace` is equivalent to a struct in C, and elements are accessed similarly (`part.density`), and can contain both
further `fspace` and other variable types. The name of the `fspace` must always be `part` to be visible to the DSL.

This is required code for all particle declarations but can easily be extended by adding additional values
into the `fspace`:
```
  fspace part{
    neighbour_part_space : neighbour_part,
    core_part_space : core_part,
    density : float,
    accel_x : float,
    accel_y : float,
    accel_z : float
  }
```

The `neighbour_part_space` is an opaque `fspace` which is used by the neighbour search algorithms, and may contain nothing. This should not be used in kernel code.

The `core_part_space` contains a set of variables that all particles are likely to have. This is declared in `src/particles/core_part.rg`,
and contains:
```
  fspace core_part{
    pos_x : double,
    pos_y : double,
    pos_z : double,
    mass : double,
    vel_x : double,
    vel_y : double,
    vel_z : double,
    cutoff : double,
    id : int1d
  }
```

These can be accessed safely in kernels through `part.core_part_space.pos_x` for example.


## Kernel Declarations
Different types of Kernel are available in HartreeParticleDSL. At the moment, pairwise and per-particle kernels are implemented.

### Pairwise Kernel

Pairwise kernels are computed on all particle pairs which are within their `cutoff` radii. Certain neighbour search algorithms may
assume a fixed global cutoff, while others will allow for per-particle cutoffs. The documentation for the neighbour search algorithms 
will discuss the specifics.

A pairwise kernel has the following declaration:
```
  function pairwise_kernel_name( part1, part2, r2 )
    local kernel = rquote
      --Kernel code goes here
    end
    return kernel
  end
```

The arguments to the function are two `part` ( `part1` and `part2` ) and `r2` which contains the square of the distance between them.
The values in the particles can be freely modified, and local variables can be created and used as required.

### Per-particle Kernel

Per-particle kernels are applied to all particles in the system. 

A per-particle kernel has the following 
declaration:
```
   function per_particle_kernel_name( part, config)
     local kernel = rquote
       --Kernel code goes here
     end
     return kernel
   end
```

The arguments to the function is a `part` and the `config` type. The `part` can be freely modified, while the `config` type is currently read-only.

### Using kernels in the main program

The DSL provides an `invoke` function, which is used to call kernels in the main program. The syntax for using the invoke functionality is:
```
[invoke(variables.config, {kernel1,type1}, {kernel2,type2}, extra options)];
```
The first argument to the invoke call is always the `variables.config` structure, defined by the DSL.

The next arguments is a group of kernels and their types. The possible types of kernels currently supported is:
1. `SYMMETRIC_PAIRWISE` - Kernels of this type are pairwise kernels that write to both particles in the interaction.
2. `ASYMMETRIC_PAIRWISE` - Kernels of this type are pairwise kernels that only write to `part1`.
3. `PER_PART` - Kernels of this type are per-particle kernels applied to all particles in the system.

Any number of kernels and types could be listed, e.g. `{kernel1, PER_PART}, ...., {kernel99, SYMMETRIC_PAIRWISE}`.

There are also a small number of extra options available for the function, which could be used by advanced users to potentially improve performance.
1. `BARRIER` - The barrier option forces a barrier to occur at the end of the invoke. This option is enabled by default.
2. `NO_BARRIER` - This removes the barrier at the end of the invoke, and will always override any `BARRIER` option chosen.
3. `MULTI_KERNEL` - This option enables the runtime system to merge multiple kernels into single tasks, which may lead to improved performance. The
runtime system will only do this if data dependencies allow, so should never cause incorrect results. This is enabled by default.
4. `SINGLE_KERNEL` - This option disables the runtime system from merging kernels into a single task. This always overrides any `MULTI_KERNEL` options.

## Main Program

The main program is broken into a few sections. 
The overall file structure would usually be similar to:
```
    import "regent"
    --DSL setup (as shown above)

    task main()
      --Code goes here
    end

    run_DSL(main)
```

This code sets up the headers and file, and the `run_DSL(main)` call starts the program on the `main` task.

Inside the `main` task there are a few section. First the code needs to initialise the data structures. This can be done manually (though not recommended), instead 
IO modules contain an initialisation function, which can be used with:
```
    [initialisation(variables, other arguments)];
```

For details on the initialisation (and finalisation or other IO functions), check the appropriate IO module's documentation.

### The timestepping loop

The main body of the method is free to be defined however you want, with the only limitation that all functions used must be either:
1. Tasks defined through the DSL's code generation
2. Explicit user-created Regent tasks
3. Code that only affects local variables

An example of this 
might be:
```
    task main()
      [initialisation(variables)];
      var time = 0.0
      var timestep = 0.001
      while(time < 1.0) do
        [invoke(...)];
        time = timestep + time
      end
    [finalisation(variables)];
    end
``` 

The final line in the file is usually `run_DSL(main)` (or `compile_DSL` to compile the code as a binary). 

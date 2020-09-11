# DSL Program Structure

RegentParticleDSL programs consist of three sections:

1. The particle type declaration
2. The kernel declaration(s)
3. The main program declaration.

Additional functionality can be added to programs, but these are the minimum required.

## The Particle Type

Regent's particle type is the main data structure used in the particle methods. A base declaration is
available in `src/particles/default_part.rg`. The outline of the particle type is:
```
  import "regent"
  require("defaults")
  
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
Different types of Kernel are available in RegentParticleDSL. At the moment, pairwise and per-particle kernels are implemented.

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

### Using kernels for code generation

Once the kernels are written, they are used with the code generation functions to create the functions that one would use in the main program.
For example, to create a per-particle function from a kernel:
```
  per_particle_function = run_per_particle_task( per_particle_kernel_name )
```

After this call, the `per_particle_function` call is usable in the main program code. For the exact functions and arguments for a specific neighbour search
algorithm, check the appropriate module's documentation.

## Main Program

The main program is broken into a few sections. 
The overall file structure would usually be similar to:
```
    import "regent"
    require("defaults")
    require("other/headers/needed")

    task main()
      --Code goes here
    end

    regentlib.start(main)
```

This code sets up the headers and file, and the `regentlib.start(main)` call starts the program on the `main` task.

Inside the `main` task there are a few section. First the code needs to initialise the data structures. At the moment this is done manually, however 
IO modules will contain an initialisation function, which can be used with:
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
    local timestepping_task = run_per_particle_task( timestep )
    local interaction_task = create_symmetric_pairwise_runner( kernel )

    task main()
      [initialisation(variables)];
      var time = 0.0
      var timestep = 0.001
      while(time < 1.0) do
        interaction_task(...)
        timestepping_task(...)
        time = timestep + time
      end
    [finalisation(variables)];
    end
``` 


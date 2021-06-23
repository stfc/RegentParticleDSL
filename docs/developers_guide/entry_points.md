# DEVELOPERS GUIDE: Entry points

The HartreeParticleDSL system provides relatively few entry points that are aimed to be used by users. This section of the documentation aims to 
describe what these are, so developers wishing to add new algorithms or replace sections of the framework should match these exactly.

## Invoke entry point

The invoke framework contains a single entry point:
```
function invoke( config, ... )
```
The `...` arguments will always contain at least one Lua table of `{function, KERNEL_TYPE}`, where the function will return a Regent `rexpr`. 
They may additionally contain some other arguments, such as whether there should be a barrier at the end of the invoke, and whether to allow
the invoke framework to combine kernels into single launches. The invoke framework should always enforce the barrier if asked, however it may
choose to ignore the user's requests for dissallowing kernel combining.

The return value from the invoke framework should be a single `rquote` which executes the code provided by the user's arguments.


## Neighbour search entry point

The neighbour search entry points are used by the user and/or by the invoke framework.

### Required system entry points.
The neighbour search algorithm must provide `neighbour_config` and `neighbour_part` field space definitions.

The `neighbour_part` may be empty, i.e.:
```
fspace neighbour_part{
}
```

The `neighbour_config` may be empty, i.e.:
```
fspace neighbour_config{
}
```

These should be invisible to the user, and users should not require knowledge of what they contain to be able to write correct code.


### User visible entry points
The neighbour search provides two entry points that are visible to the user to initialise the neighbour search algorithm:
1. `neighbour_init.initialise(variables)`
2. `neighbour_init.update_cells(variables)`

Both of these functions should return `rquote` which contain any initialisation code required for the neighbour search algorithm to run correctly.

The functions are allowed to overwrite what `variables` fields point to, and are single-entry points for a given `variables`.

### Invoke visible entry points

The neighbour search should provide 4 entry points to the invoke system:
1. `create_symmetric_pairwise_runner( kernel_name, config )`
2. `create_asymmetric_pairwise_runner( kernel_name, config )`
3. `run_per_particle_task( kernel_name, config )`

and 

```
__demand(__inline)
task neighbour_init.check_valid(ns : neighbour_part)
```

The first 3 should generate code to run the kernels provided for all of the particles in the system, and return an `rquote` containing all of the required code. These should be re-entrant and may be called many times throughout the program.

The final one should return a boolean to tell a caller whether a given particle is a real particle or a dummy particle used by the neighbour search algorithm to avoid
repartitioning etc. If no dummy particles exist, then this can just be:
```
__demand(__inline)
task neighbour_init.check_valid(ns: neighbour_part)

    return true
end
```

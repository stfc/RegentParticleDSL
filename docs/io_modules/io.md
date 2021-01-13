# Io Modules

Only one IO module has been implemented so far. IO modules should not be included in the `default.rg`, and should be the
last file in the `require` list in the main program file.

## Empty IO module

The empty IO module doesn't use File IO, but provides an initialisation function to enable the creation
of particles manually. The module is imported with:
```
require("src/io_modules/empty_io/import_empty_io")
```


 The initialisation function has the following signature:
```
[initialisation_function(variables, part_count, space_dimension_x, space_dimension_y, space_dimension_z)];
```
Where `variables` is an included set of variable names from the IO module, the `part_count` as a user specified
number of particles in the simulation, and the `space_dimension_{x,y,z}` values contain the size of the 
simulation domain (which at current must be cuboidal).

The Empty IO module provides no finalisation function as the DSL does not require cleanup without file IO accesses.

### Creating a particle list
After using the initialisation function, all of the particles `core_part_space` values will be set to 0, and any user-defined
particle values will be uninitialised. To initialise them, there needs to be a user-defined task to set value. For example:
```
task init_particles( particle_array : region(ispace(int1d), part)) where writes(particle_array) do
  --Create a counter to set the particle IDs
  var counter : int = 0
  --For this example code, we know we have 125 particles, and we want a 5x5x5 grid
  for x= 0,5 do
    for y = 0,5 do
      for z = 0,5 do
        particle_array[counter].core_part_space.pos_x = double(x) * 0.2
        particle_array[counter].core_part_space.pos_y = double(y) * 0.2
        particle_array[counter].core_part_space.pos_z = double(z) * 0.2
        particle_array[counter].core_part_space.id = int1d(counter)
        counter = counter + 1
      end
    end
  end
end
``` 

## Simple HDF5 Module
The simple HDF5 module enables HDF5 file IO through a user-defined mapping between the HDF5 paths and the particle type.

The module is imported with
```
simple_hdf5_module = require("src/io_modules/HDF5/HDF5_simple_module")
```
and the functions are all encapsulated inside `simple_hdf5_module`.

The include path of the HDF5 library needs to be set in the `HDF5_INCLUDE_PATH` environment variable. 
If using the RegentParticleDSL docker container this will already be set.

### File Input
The simple HDF5 module provides an initialisation function to start the program from a HDF5 file:
```
simple_hdf5_module.initialisation( filename, mapper, variables, space_x, space_y, [space_z])
```
Parameters:
- `filename`: A full or relative path to the HDF5 input file.
- `mapper`: A HDF5-compatible mapper between the HDF5 file and the particle type. This is discussed later in this documentation.
- `variables`: The variables table containing the symbols used in the application.
- `space_x`: The size of the x dimension of the simulation space.
- `space_y`: The size of the y dimension of the simulation space.
- `space_z`: _optional_. The size of the z dimension of the simulation space. If not specified the simulation is assumed to be 2 dimensional.

Usage:
The `initialisation` function is called from the `main` task in the user defined program,
and requires surrounding square brackets([...];). For example:
```
task main()
  [simple_hdf5_module.initialisation(file, mapper, variables.particle_array, 1.0, 1.0, 1.0);
...
end
```
would call the initialisation function using the path declared in the `file` variable, using the mapper defined in the `mapper` variable.`


### Setting up File Ouput
Before files can be written by the HDF5 module, this function must be called:
```
simple_hdf5_module.initialise_io_module( particle_array, mapper )
```
Parameters:
- `particle_array`: The initial particle array symbol. Usually this is `variables.particle_array`.
- `mapper` : A HDF5-compatible mapper between the HDF5 file format and the particle type. This is discussed later in this documentation.

This function should be called exactly once (so before the main loop) before any file output calls, and after the particle array is initialised
(so after the `initialisation` call or similar).

### File Output
The simple HDF5 module provides a function to output the current simulation state to a HDF5 file:
```
simple_hdf5_module.write_output( filename, mapper, particle_array )
```

Parameters:
- `filename`: A full or relative path to the required HDF5 output file location. This must be a regentlib.string variable (discussed shortly).
- `mapper`: A HDF5-compatible mapper between the HDF5 file and the particle type. This is discussed later in this documentation.
- `particle_array`: The particle array symbol. Usually this is `variables.particle_array`


Usage:
The `write_output` function is called from the `main` task in the user defined program,
and requires surrounding square brackets([...];). For example:
```
task main()
...
[simple_hdf5_module.write_output( filename, mapper, variables.particle_array)];
```
would create a new hdf5 file at the location stored in `filename` using the mapper defined in the `mapper` variable.

The filename variable must be a `regentlib.string` typed variable. The easiest way to do this (for a fixed string) is to simply declare it:
```
var filename : regentlib.string = [regentlib.string]("path/to/file.hdf5")
```

If you want a filepath that varies (for example with step count), you can instead do:
```
var rawfile : rawstring = [rawstring] regentlib.c.malloc(1024)
format.snprint(rawfile, 1024, "path/to/file_{}.hdf5", step);
var filename : regentlib.string = [regentlib.string](rawfile)
[simple_hdf5_module.write_output(filename, ...)];
regentlib.c.free(rawfile)
```



### Mappers for the simple HDF5 module
The mapper is a required part of the HDF5 module, which enables the program to understand how to read and write HDF5 files for a given application.

The mapper usually needs to be declared in the lua section of the user code
 (i.e. outside of the main task), 
though other IO modules may have their own definitions to enable conversion between file formats.

The mapper is a lua table, and contains a map from "HDF5Name" to "part_field_name", for example:
```
test_mapper = {} --Define the mapper
test_mapper["Position_x"] = "core_part_space.pos_x"
test_mapper["Position_y"] = "core_part_space.pos_y"
test_mapper["Position_z"] = "core_part_space.pos_z"
```

This example maps the HDF5 paths `Position_{x,y,z}` to the `particle_array.core_part_space.{pos_x,y,z}`
sections of the particle structure.

NOTE: At current the simple HDF5 module does not support more than one nested field space,
 see [this issue](https://github.com/stfc/RegentParticleDSL/issues/41) and let us know
 if this is required functionality so we can add it.

### HDF5 read file
The simple HDF5 module also provides a function to read a file to an already existing region. This is primarily used for testing purposes, but advanced users may also
want this available.
```
simple_hdf5_module.read_file(filename, mapper, particle_array)
```

Parameters:
- `filename`: A full or relative path to the HDF5 input file.
- `mapper`: A HDF5-compatible mapper between the HDF5 file and the particle type.
- `particle_array`: The particle array symbol.

Usage:
The `write_output` function is called from the `main` task in the user defined program,
and requires surrounding square brackets([...];). For example:
```
task main()
...
[simple_hdf5_module.read_file( "path/to/input.hdf5", mapper, variables.particle_array2)];
```
would read the hdf5 file at `path/to/input.hdf5` using the mapper defined in the `mapper` variable.

## ISPH IO Module
The ISPH IO module handles the ISPH file format, and provides some other helpful functionality for ISPH files.

The module has a 2D and 3D version, and they can be imported with
```
2d_isph_module = require("src/io_modules/ISPH/2d_isph_module")
3d_isph_module = require("src/io_modules/ISPH/3d_isph_module")
```
and the functions are all encapsulated inside `XX_isph_module`.

### Initialisation
The ISPH module provides an initialisation function, which initalises a particle array from and ISPH formatted file. For the 2D version this is:
```
isph_module.initialisation_function(filename, variables, space_x, space_y, [mapper])
``` 
Parameters:
- `filename`: A full or relative path to the ISPH input file.
- `variables`: The variables definition.
- `space_x`: The size of the x dimension of the simulation space.
- `space_y`: The size of the y dimension of the simulation space.
- `mapper`: _optional_ An ISPH mapper. If not defined the function uses the default definition, which requires the particle to contain a `core_part`, and two `double` fields named `pressure` and `volume`.

The 3D version this is:
```
isph_module.initialisation_function(filename, variables, space_x, space_y, space_z, [mapper])
```
Parameters:
- `filename`: A full or relative path to the ISPH input file.
- `variables`: The variables definition.
- `space_x`: The size of the x dimension of the simulation space.
- `space_y`: The size of the y dimension of the simulation space.
- `space_z`: The size of the z dimension of the simulation space.
- `mapper`: _optional_ An ISPH mapper. If not defined the function uses the default definition, which requires the particle to contain a `core_part`, and two `double` fields named `pressure` and `volume`.

Usage:
The `initialisation` function is called from the `main` task in the user defined program,
and requires surrounding square brackets([...];). For example:
```
task main()
...
[isph_module.initialisation_function( "path/to/input.txt", variables, 3.0, 3.0)];
```
would initialise the application using the `path/to/input.txt` file, and create a space of dimensions 3.0 x 3.0.

### ISPH file output
The ISPH module also provides the ability to write output similar to the ISPH format:
```
isph_module.write_output(filename, variables, [mapper])
```

Parameters:
- `filename`: A full or relative path to the ISPH input file.
- `variables`: The variables definition.
- `mapper`: _optional_ An ISPH mapper. If not defined the function uses the default definition, which requires the particle to contain a `core_part`, and two `double` fields named `pressure` and `volume`.


Usage:
The `initialisation` function is called from the `main` task in the user defined program,
and requires surrounding square brackets([...];). For example:
```
task main()
...
[isph_module.write_output( "path/to/output.txt", variables)];
```
would create a new ISPH file at `path/to/output.txt` and store the values from the simulation.

### The ISPH Mapper

The mapper is an optional part of the ISPH module, which enables the user to specify how to store the values read from an ISPH file.

If used, the mapper needs to be declared in the lua section of the user code
 (i.e. outside of the main task).

The mapper is a lua table, and contains a map from the column in the file to the part structure field, for example, the inbuilt default mapper is defined for the 2D case as:
```
local isph_mapper = {}
isph_mapper[1] = "core_part_space.pos_x"
isph_mapper[2] = "core_part_space.pos_y"
isph_mapper[3] = "core_part_space.vel_x"
isph_mapper[4] = "core_part_space.vel_y"
isph_mapper[5] = "pressure"
isph_mapper[6] = "volume"
```

and for the 3D case as:
```
local isph_mapper = {}
isph_mapper[1] = "core_part_space.pos_x"
isph_mapper[2] = "core_part_space.pos_y"
isph_mapper[3] = "core_part_space.pos_z"
isph_mapper[4] = "core_part_space.vel_x"
isph_mapper[5] = "core_part_space.vel_y"
isph_mapper[6] = "core_part_space.vel_z"
isph_mapper[7] = "pressure"
isph_mapper[8] = "volume"
```

This mapper stores the first value in each particle's definition from the ISPH file into `core_part_space.pos_x`, the second in `core_part_space.pos_y` and so on.

### ISPH to HDF5 conversion
In addition to the ISPH file functionality, the ISPH module contains a mapper to enable outputting HDF5 files using the simple HDF5 module. This mapper is defined as `isph_module.hdf5_mapper`, and can be passed to the simple HDF5 module functions to enable HDF5 IO with ISPH file formats.

For example, this could be used to create a simple conversion between the file formats:
```
task main()
[isph_module.initialisation_function("/path/to/input.txt", variables, 1.0, 1.0)];
[hdf5_module.write_output("/path/to/output.hdf5", isph_module.hdf5_mapper, variables.particle_array)];
end
```

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

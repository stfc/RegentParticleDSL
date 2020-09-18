import "regent"

require("defaults")
require("src/particles/init_part")
require("src/neighbour_search/cell_pair/neighbour_search")

format = require("std/format")
local c = regentlib.c

isph_module = {}

local stdio = terralib.includec("stdio.h")
local stdlib = terralib.includec("stdlib.h")

--Hack to deal with the fact that the EOF defition in stdio.h is not imported
--TODO: Improve this code because its linux-dependent right now...
local EOF = -1

local task read_file(filename : rawstring, particle_region : region(ispace(int1d), part)) where writes(particle_region) do
  var f = stdio.fopen(filename, "r")
  for part in particle_region.ispace do
    var pos_x : double
    var vel_x : double
    var pos_y : double
    var vel_y : double
    var rho : double
    var who_knows : double
    stdio.fscanf(f, "%lf %lf %lf\n %lf %lf %lf\n", &vel_x, &pos_x, &vel_y, &pos_y, &rho, &who_knows)
    particle_region[part].core_part_space.pos_x = pos_x
    particle_region[part].core_part_space.vel_x = vel_x
    particle_region[part].core_part_space.pos_y = pos_y
    particle_region[part].core_part_space.vel_y = vel_y
    --TODO Make use of the remaining data
  end
  stdio.fclose(f)
end

local task write_file(filename : rawstring, particle_region : region(ispace(int1d), part)) where reads(particle_region) do
  var f = stdio.fopen(filename, "w")
  for part in particle_region.ispace do
    stdio.fprintf(f, "%e %e %e\n %e %e %e\n", particle_region[part].core_part_space.vel_x,
                                           particle_region[part].core_part_space.pos_x,
                                           particle_region[part].core_part_space.vel_y,
                                           particle_region[part].core_part_space.pos_y,
                                           particle_region[part].core_part_space.vel_z,
                                           particle_region[part].core_part_space.pos_z)
               
  end
  stdio.fclose(f)
end

local task particle_initialisation(particle_region : region(ispace(int1d), part), filename : rawstring) where writes(particle_region) do

zero_core_part(particle_region)
zero_neighbour_part(particle_region)
read_file(filename, particle_region)
end

local terra read_part_count(filename : rawstring) : uint32
  var f = stdio.fopen(filename, "r")
  var lines : uint32 = 0
  while(EOF ~= ( stdio.fscanf(f, "%*[^\n]"))) do
    stdio.fscanf(f, "%*c")
    lines = lines + 1
  end
  stdio.fclose(f)
  return lines
end

function isph_module.initialisation_function(filename, variables, space_x, space_y)

local init_string = rquote
  regentlib.assert(space_x > 0.0, "Space x dimension must be greater than 0")
  regentlib.assert(space_y > 0.0, "Space y dimension must be greater than 0")
  var space_z = 0.0
  var part_count = read_part_count(filename) / 2
  format.println("part_count {}", part_count)
  var particles_space = ispace(int1d, part_count)
  var [variables.particle_array] = region(particles_space, part)
  var [variables.config] = region(ispace(int1d, 1), config_type)
  fill([variables.config].{space.dim_x, space.dim_y, space.dim_z}, 0.0)

  init_space(space_x, space_y, space_z, [variables.config])
  particle_initialisation([variables.particle_array], filename)
  c.legion_runtime_issue_execution_fence(__runtime(), __context())
end
return init_string
end

function isph_module.write_output(filename, variables)
local write_string = rquote
write_file(filename, [variables.particle_array])
end
return write_string
end

--fspace hdf5_io_space{
--  Position_x : double,
--  Position_y : double,
--  Velocity_x : double,
--  Velocity_y : double
--}

isph_module.hdf5_mapper = {}
isph_module.hdf5_mapper["Position_x"] = "core_part_space.pos_x"
isph_module.hdf5_mapper["Position_y"] = "core_part_space.pos_y"
isph_module.hdf5_mapper["Velocity_x"] = "core_part_space.vel_x"
isph_module.hdf5_mapper["Velocity_y"] = "core_part_space.vel_y"

return isph_module

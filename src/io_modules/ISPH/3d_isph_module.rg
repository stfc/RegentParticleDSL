-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

require("src/particles/init_part")

string_to_field_path = require("src/utils/string_to_fieldpath")
format = require("std/format")
local c = regentlib.c

isph_module = {}

local stdio = terralib.includec("stdio.h")
local stdlib = terralib.includec("stdlib.h")

--Hack to deal with the fact that the EOF defition in stdio.h is not imported
--TODO: Improve this code because its linux-dependent right now...
--https://github.com/stfc/RegentParticleDSL/issues/39
local EOF = -1

--The default isph mapper - ISPH files always have 6 values to read in, and this is used
--if the user doesn't define their own mapper.
local isph_mapper = {}
isph_mapper[1] = "core_part_space.pos_x"
isph_mapper[2] = "core_part_space.pos_y"
isph_mapper[3] = "core_part_space.pos_z"
isph_mapper[4] = "core_part_space.vel_x"
isph_mapper[5] = "core_part_space.vel_y"
isph_mapper[5] = "core_part_space.vel_z"
isph_mapper[7] = "pressure"
isph_mapper[8] = "volume"

--Generate the read_file_task from the mapper
local function gen_read_file_task(mapper)
  --Create the field paths from the mapper
  local field_1 = string_to_field_path.get_field_path( mapper[1] ) 
  local field_2 = string_to_field_path.get_field_path( mapper[2] ) 
  local field_3 = string_to_field_path.get_field_path( mapper[3] ) 
  local field_4 = string_to_field_path.get_field_path( mapper[4] ) 
  local field_5 = string_to_field_path.get_field_path( mapper[5] ) 
  local field_6 = string_to_field_path.get_field_path( mapper[6] ) 
  local field_7 = string_to_field_path.get_field_path( mapper[7] ) 
  local field_8 = string_to_field_path.get_field_path( mapper[8] ) 
  local task read_file(filename : rawstring, particle_region : region(ispace(int1d), part)) where writes(particle_region) do
    var f = stdio.fopen(filename, "r")
    for part in particle_region.ispace do
      var val_1 : double
      var val_2 : double
      var val_3 : double
      var val_4 : double
      var val_5 : double
      var val_6 : double
      var val_7 : double
      var val_8 : double
      --Read in 6 values across 2 lines, and store them in the mapper defined fields in the particle structures
      stdio.fscanf(f, "%lf %lf %lf\n %lf %lf %lf\n %lf %lf\n", &val_1, &val_2, &val_3, &val_4, &val_5, &val_6, &val_7, &val_8)
      particle_region[part].[field_1] = val_1
      particle_region[part].[field_2] = val_2
      particle_region[part].[field_3] = val_3
      particle_region[part].[field_4] = val_4
      particle_region[part].[field_5] = val_5
      particle_region[part].[field_6] = val_6
      particle_region[part].[field_7] = val_7
      particle_region[part].[field_8] = val_8
    end
    stdio.fclose(f)
  end
  return read_file
end

--Generate the write_file_task from the mapper
local function gen_write_file_task(mapper)
  --Create the field paths from the mapper
  local field_1 = string_to_field_path.get_field_path( mapper[1] ) 
  local field_2 = string_to_field_path.get_field_path( mapper[2] ) 
  local field_3 = string_to_field_path.get_field_path( mapper[3] ) 
  local field_4 = string_to_field_path.get_field_path( mapper[4] ) 
  local field_5 = string_to_field_path.get_field_path( mapper[5] ) 
  local field_6 = string_to_field_path.get_field_path( mapper[6] ) 
  local field_6 = string_to_field_path.get_field_path( mapper[7] ) 
  local field_6 = string_to_field_path.get_field_path( mapper[8] ) 
  local task write_file(filename : rawstring, particle_region : region(ispace(int1d), part)) where reads(particle_region) do
    var f = stdio.fopen(filename, "w")
    for part in particle_region.ispace do
      --Write out 6 values in the ISPH format, according to the mapper defined fields
      stdio.fprintf(f, "%e %e %e\n %e %e %e\n %e %e\n", particle_region[part].[field_1],
                                             particle_region[part].[field_2],
                                             particle_region[part].[field_3],
                                             particle_region[part].[field_4],
                                             particle_region[part].[field_5],
                                             particle_region[part].[field_6],
                                             particle_region[part].[field_7],
                                             particle_region[part].[field_8])
                 
    end
    stdio.fclose(f)
  end
  return write_file
end


local particle_initialisation = generate_zero_part_func()

--Read in the number of particles. This is the number of lines in the file / 2
local terra read_part_count(filename : rawstring) : uint32
  var f = stdio.fopen(filename, "r")
  var lines : uint32 = 0
  while(EOF ~= ( stdio.fscanf(f, "%*[^\n]"))) do
    stdio.fscanf(f, "%*c")
    lines = lines + 1
  end
  stdio.fclose(f)
  return lines / 2
end

--The isph module's initialisation function. This is called from user code, and
--initialises the particle structures according to the main library code and the user's mapper.
--Mapper argument is optional, if not given then the default isph mapper is used
function isph_module.initialisation_function(filename, variables, space_x, space_y, space_z, mapper)
  if mapper == nil then
    mapper = isph_mapper
  end
  local read_task = gen_read_file_task(mapper)
local init_string = rquote
  regentlib.assert(space_x > 0.0, "Space x dimension must be greater than 0")
  regentlib.assert(space_y > 0.0, "Space y dimension must be greater than 0")
  regentlib.assert(space_z > 0.0, "Space y dimension must be greater than 0")
  var part_count = read_part_count(filename)
  format.println("Reading {} ISPH particles", part_count)
  var particles_space = ispace(int1d, part_count)
  var [variables.particle_array] = region(particles_space, part)
  var [variables.config] = region(ispace(int1d, 1), config_type)
  fill([variables.config].{space.dim_x, space.dim_y, space.dim_z}, 0.0)

  init_space(space_x, space_y, space_z, [variables.config])
--  particle_initialisation([variables.particle_array], filename)
  particle_initialisation([variables.particle_array])
  read_task(filename, [variables.particle_array])
  c.legion_runtime_issue_execution_fence(__runtime(), __context())
end
return init_string
end

--The isph module's file writing function. This is called from user code and
--writes filename according to the mapper.
--Mapper argument is optional, if not given then the default isph mapper is used
function isph_module.write_output(filename, variables, mapper)
  if mapper == nil then
    mapper = isph_mapper
  end
  local write_task = gen_write_file_task(mapper)

local write_string = rquote
write_task(filename, [variables.particle_array])
end
return write_string
end

--Default ISPH to HDF5 mapper, can be used by the user to create HDF5
--output files from ISPH particle systems.
isph_module.hdf5_mapper = {}
isph_module.hdf5_mapper["Position_x"] = "core_part_space.pos_x"
isph_module.hdf5_mapper["Position_y"] = "core_part_space.pos_y"
isph_module.hdf5_mapper["Position_z"] = "core_part_space.pos_z"
isph_module.hdf5_mapper["Velocity_x"] = "core_part_space.vel_x"
isph_module.hdf5_mapper["Velocity_y"] = "core_part_space.vel_y"
isph_module.hdf5_mapper["Velocity_z"] = "core_part_space.vel_z"
isph_module.hdf5_mapper["Pressure"] = "pressure"
isph_module.hdf5_mapper["Volume"] = "volume"

return isph_module

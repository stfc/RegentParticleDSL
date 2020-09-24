-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

require("defaults")
local format = require("std/format")
--TODO: We want to make this not just specific to a single Issue: #46
require("src/particles/core_part")
require("src/neighbour_search/cell_pair_v2/import_cell_pair")

require("defaults")
require("src/neighbour_search/cell_pair_v2/import_cell_pair")
require("src/neighbour_search/cell_pair_v2/neighbour_search")
require("src/neighbour_search/cell_pair_v2/cell")
require("examples/interaction_count/interaction_count_kernel")
require("examples/interaction_count/infrastructure/interaction_count_init")
simple_hdf5_module = require("src/io_modules/HDF5/HDF5_simple_module")

variables = {}

variables.config = regentlib.newsymbol("config")
variables.particle_array = regentlib.newsymbol("particle_region")
variables.solution_array = regentlib.newsymbol("solution_region")
variables.cell_space = regentlib.newsymbol("cell_space")

local c = regentlib.c
format = require("std/format")
local stdlib = terralib.includec("stdlib.h")

local hdf5_read_mapper = {}
hdf5_read_mapper["cutoff"] = "core_part_space.cutoff"
hdf5_read_mapper["pos_x"] = "core_part_space.pos_x"
hdf5_read_mapper["pos_y"] = "core_part_space.pos_y"
hdf5_read_mapper["pos_z"] = "core_part_space.pos_z"
hdf5_read_mapper["ids"] = "core_part_space.id"

local hdf5_write_mapper = {}
hdf5_write_mapper["cutoff"] = "core_part_space.cutoff"
hdf5_write_mapper["interactions"] = "interactions"
hdf5_write_mapper["pos_x"] = "core_part_space.pos_x"
hdf5_write_mapper["pos_y"] = "core_part_space.pos_y"
hdf5_write_mapper["pos_z"] = "core_part_space.pos_z"
hdf5_write_mapper["ids"] = "core_part_space.id"

local c = regentlib.c
local cstring = terralib.includec("string.h")

local input_file = "tests/interaction_count/test.hdf5"
local solution_file = "tests/interaction_count/test2.hdf5"
local x_cell = 3.0
local y_cell = 3.0
local z_cell = 3.0


local function getarg(key)
  for i=0, #arg-1 do
    if(arg[i] == key) then
      return arg[i+1]
    end
  end
  print("FAILURE: Failed to find input argument")
  os.exit(1)
  return nil
end

local function get_optional_arg(key)
  for i=0, #arg-1 do
    if(arg[i] == key) then
      return arg[i+1]
    end
  end
  return nil
end

local function read_args()
  local read_val = get_optional_arg("-input")
  if(read_val ~= nil) then
    input_file = read_val
  end
  read_val = get_optional_arg("-solution")
  if(read_val ~= nil) then
    solution_file = read_val
  end
  read_val = get_optional_arg("-x_cell")
  if(read_val ~= nil) then
    x_cell = tonumber(read_val)
  end
  read_val = get_optional_arg("-y_cell")
  if(read_val ~= nil) then
    y_cell = tonumber(read_val)
  end
  read_val = get_optional_arg("-z_cell")
  if(read_val ~= nil) then
    z_cell = tonumber(read_val)
  end
end
read_args()



--Symmetric interaction count kernel
function symmetric_interaction_count_kernel(part1, part2, r2)
local kernel = rquote
  if part1.core_part_space.id == int1d(0) or part2.core_part_space.id == int1d(0) then
    format.println("{} and {}", part1.core_part_space.id, part2.core_part_space.id)
  end
  part1.interactions = part1.interactions + 1
  part2.interactions = part2.interactions + 1
end
return kernel
end

local interaction_tasks_runner = create_symmetric_pairwise_runner( symmetric_interaction_count_kernel, variables.config, variables.cell_space )

task comparison(computed : region(ispace(int1d), part), solution : region(ispace(int1d), part)) where
reads(computed.interactions, solution.interactions, computed.core_part_space.id, solution.core_part_space.id) do
  for x in computed.ispace do
    for y in solution.ispace do
      if computed[x].core_part_space.id == solution[x].core_part_space.id then
        if computed[x].interactions ~= solution[x].interactions then
          var buffer = [rawstring](regentlib.c.malloc(1024))
          format.snprintln(buffer, 1024, "Interactions incorrect for particle ID {}, computed {} solution {}", computed[x].core_part_space.id,
                            computed[x].interactions, solution[x].interactions)
          regentlib.assert(computed[x].interactions == solution[x].interactions, buffer)
          regentlib.c.free(buffer)
        end
      end
    end
  end
  format.println("All interactions computed correctly!")
end

task main_task()
  [simple_hdf5_module.initialisation( input_file, hdf5_read_mapper, variables, x_cell, y_cell, z_cell)];
--
  particles_to_cell_launcher(variables.particle_array, variables.config);
  var [variables.cell_space] = update_cell_partitions(variables.particle_array, variables.config);
   [interaction_tasks_runner];

  [simple_hdf5_module.read_file( solution_file, hdf5_write_mapper, variables.solution_array)];
  comparison(variables.particle_array, variables.solution_array);
end

regentlib.start(main_task)

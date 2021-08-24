-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

require("src/RegentParticleDSL")
set_dimensionality(3)
set_periodicity(true)
--disable_high_performance()
setup_part()
local format = require("std/format")
--TODO: We want to make this not just specific to a single Issue: #46
require("src/particles/core_part")
require("examples/interaction_count/interaction_count_part")
setup_dsl()

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

--local input_file = "tests/interaction_count/test.hdf5"
--local solution_file = "tests/interaction_count/test2.hdf5"
--local x_cell = 3.0
--local y_cell = 3.0
--local z_cell = 3.0
--
--
--local function getarg(key)
--  for i=0, #arg-1 do
--    if(arg[i] == key) then
--      return arg[i+1]
--    end
--  end
--  print("FAILURE: Failed to find input argument")
--  os.exit(1)
--  return nil
--end
--
--local function get_optional_arg(key)
--  for i=0, #arg-1 do
--    if(arg[i] == key) then
--      return arg[i+1]
--    end
--  end
--  return nil
--end
--
--local function read_args()
--  local read_val = get_optional_arg("-input")
--  if(read_val ~= nil) then
--    input_file = read_val
--  end
--  read_val = get_optional_arg("-solution")
--  if(read_val ~= nil) then
--    solution_file = read_val
--  end
--  read_val = get_optional_arg("-x_cell")
--  if(read_val ~= nil) then
--    x_cell = tonumber(read_val)
--  end
--  read_val = get_optional_arg("-y_cell")
--  if(read_val ~= nil) then
--    y_cell = tonumber(read_val)
--  end
--  read_val = get_optional_arg("-z_cell")
--  if(read_val ~= nil) then
--    z_cell = tonumber(read_val)
--  end
--end
--read_args()

local cstring = terralib.includec("string.h")
local cstdio = terralib.includec("stdio.h")
local terra get_optional_arg( arg_name : rawstring )
  var args = c.legion_runtime_get_input_args()
  var i = 1
  while i < args.argc do
--    cstdio.printf("%s\n", args.argv[i])
    if cstring.strcmp(args.argv[i], arg_name) == 0 then
      if i + 1 < args.argc then
        return args.argv[i+1]
      else
        return nil
      end
    end
    i = i + 1
  end
  return nil
end

--Aymmetric interaction count kernel
--function asymmetric_interaction_count_kernel(part1, part2, r2)
--local kernel = rquote
--  part1.interactions = part1.interactions + 1
--end
--return kernel
--end
function asymmetric_interaction_count_kernel(part1, part2, r2)
local kernel = rquote
  part1.interactions += 1
end
return kernel
end

function move_part(part)
local kernel = rquote
    part.core_part_space.pos_x += 0.25
    part.core_part_space.pos_y += 0.15
    part.core_part_space.pos_z += 0.05
    part.interactions = 0
end
return kernel
end

task comparison(computed : region(ispace(int1d), part), solution : region(ispace(int1d), part), step : int) where
reads(computed.interactions, solution.interactions, computed.core_part_space.id, solution.core_part_space.id, computed.neighbour_part_space) do
  for x in computed.ispace do
    for y in solution.ispace do
      if computed[x].core_part_space.id == solution[y].core_part_space.id and computed[x].neighbour_part_space._valid then
        if computed[x].interactions ~= solution[y].interactions then
          var buffer = [rawstring](regentlib.c.malloc(1024))
          format.snprintln(buffer, 1024, "Step {}: Interactions incorrect for particle ID {}, computed {} solution {}", step, computed[x].core_part_space.id,
                            computed[x].interactions, solution[x].interactions)
          regentlib.assert(computed[x].interactions == solution[x].interactions, buffer)
          regentlib.c.free(buffer)
        end
      end
    end
  end
  format.println("All interactions computed correctly!")
end

local input_file = regentlib.newsymbol("input_file")
local solution_file = regentlib.newsymbol("solution_file")
local x_cell = regentlib.newsymbol("x_cell")
local y_cell = regentlib.newsymbol("y_cell")
local z_cell = regentlib.newsymbol("z_cell")

task main_task()
  var [input_file] = [regentlib.string](get_optional_arg("-input"));
  var [solution_file] = [regentlib.string](get_optional_arg("-solution"));
  var [x_cell] = c.atof(get_optional_arg("-x_cell"));
  var [y_cell] = c.atof(get_optional_arg("-y_cell"));
  var [z_cell] = c.atof(get_optional_arg("-z_cell"));
  [simple_hdf5_module.initialisation( input_file, hdf5_read_mapper, variables, x_cell, y_cell, z_cell)];
--
  [neighbour_init.initialise(variables)];
  [neighbour_init.update_cells(variables)];
  [simple_hdf5_module.read_file( solution_file, hdf5_write_mapper, variables.solution_array)];


  for i=0, 100 do
    [invoke(variables.config, {asymmetric_interaction_count_kernel,ASYMMETRIC_PAIRWISE}, BARRIER)]; 
    comparison(neighbour_init.padded_particle_array, variables.solution_array, i);
    [invoke(variables.config, {move_part, PER_PART}, BARRIER)];

  end
end

terra set_mappers()

end

  compile_DSL( main_task, "tests/moving_interaction_count/interaction_test_asym.exe")

  --local root_dir = "./tests/interaction_count/"
  --local out_dir = (os.getenv('OBJNAME') and os.getenv('OBJNAME'):match('.*/')) or root_dir
  --local link_flags = terralib.newlist({"-L" .. out_dir, "-lm", "-lhdf5"})
  --local exe = os.getenv('OBJNAME') or "tests/interaction_count/interaction_test_asym.exe"
  --regentlib.saveobj(main_task, exe, "executable", set_mappers, link_flags)
--regentlib.start(main_task)
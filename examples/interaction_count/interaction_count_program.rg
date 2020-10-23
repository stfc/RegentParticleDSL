-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

require("defaults")
require("src/neighbour_search/cell_pair_tradequeues/import_cell_pair")
neighbour_init = require("src/neighbour_search/cell_pair_tradequeues/neighbour_init")
require("src/neighbour_search/cell_pair_tradequeues/neighbour_search")
require("src/neighbour_search/cell_pair_tradequeues/cell")
require("examples/interaction_count/interaction_count_kernel")
require("examples/interaction_count/infrastructure/interaction_count_init")
simple_hdf5_module = require("src/io_modules/HDF5/HDF5_simple_module")

local c = regentlib.c
format = require("std/format")
variables = require("examples/interaction_count/infrastructure/interaction_count_variables")
local stdlib = terralib.includec("stdlib.h")

local hdf5_read_mapper = {}
hdf5_read_mapper["cutoff"] = "core_part_space.cutoff"
hdf5_read_mapper["pos_x"] = "core_part_space.pos_x"
hdf5_read_mapper["pos_y"] = "core_part_space.pos_y"
hdf5_read_mapper["pos_z"] = "core_part_space.pos_z"

local hdf5_write_mapper = {}
hdf5_write_mapper["cutoff"] = "core_part_space.cutoff"
hdf5_write_mapper["interactions"] = "interactions"
hdf5_write_mapper["pos_x"] = "core_part_space.pos_x"
hdf5_write_mapper["pos_y"] = "core_part_space.pos_y"
hdf5_write_mapper["pos_z"] = "core_part_space.pos_z"
hdf5_write_mapper["ids"] = "core_part_space.id"


task main_task()

  --We will use 9 particles, with a large cutoff. The initialisation is currently declared in the interaction_count_init file.
  --In the future much of this will be abstracted into the DSL, though user-specified particles can be done easily enough, and
  --example code to do this will be added.
  [simple_hdf5_module.initialisation("examples/interaction_count/basic_test.hdf5", hdf5_read_mapper, variables, 3.0, 3.0, 3.0)];
  for x in [variables.particle_array].ispace do
    [variables.particle_array][x].core_part_space.id = int1d(x)
  end
  
   [neighbour_init.initialise(variables)];
  [neighbour_init.update_cells(variables)];
 
  --Run the interaction tasks "timestep". We could do this multiple types and know it will be correct due to 
  --the programming model
   [invoke(variables.config, {symmetric_interaction_count_kernel, SYMMETRIC_PAIRWISE}, BARRIER)];
--   [invoke(variables.config, {asymmetric_interaction_count_kernel, SYMMETRIC_PAIRWISE}, BARRIER)];

  --This finalisation function is declared in the interaction_count_init file.
  --It contains various tests for the basic test to check correctness.
  [finalisation_function(neighbour_init.padded_particle_array, variables.config)];
  [simple_hdf5_module.write_output("examples/interaction_count/basic_test.hdf5", hdf5_write_mapper, neighbour_init.padded_particle_array)];
end

regentlib.start(main_task)

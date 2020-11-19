-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

require("defaults")
require("src/neighbour_search/cell_pair_tradequeues/cell")
require("src/particles/init_part")

local c = regentlib.c

--TODO: This is hopefully what the initialisation function looks like for now
--It can of course read from files instead when necessary.
function initialisation_function(particle_array, config) 

local init_string = rquote
  var particles_space = ispace(int1d, 9)
  var [particle_array] = region(particles_space, part)
  var [config] = region(ispace(int1d, 1), config_type)
  fill([config].{space.dim_x, space.dim_y, space.dim_z}, 0.0)
  
  
  init_space(3.0, 3.0, 3.0, [config])
  particle_initialisation([particle_array])
  initialise_cells([config], [particle_array])
end
return init_string
end

--This currently is wrong
function finalisation_function(particle_array, config)

local final_string = rquote
  c.legion_runtime_issue_execution_fence(__runtime(), __context())
  format.println("{}", [particle_array][0].interactions)
  for point in [particle_array].ispace do
    if [particle_array][point].neighbour_part_space._valid then
      var s : rawstring
      s = [rawstring] (regentlib.c.malloc(256))
      format.snprintln(s,256, "test failed for {} with id {}: value {}", point, [variables.particle_array][point].core_part_space.id, [variables.particle_array][point].interactions)
      regentlib.assert([variables.particle_array][point].interactions == 8, s)
      regentlib.c.free(s)
    end
  end
  
 -- write_hdf5_snapshot("examples/interaction_count/basic_test.hdf5", [particle_array], [config])
  c.legion_runtime_issue_execution_fence(__runtime(), __context())
end
return final_string
end

local zero_part_task = generate_zero_part_func()
task particle_initialisation(particle_region : region(ispace(int1d), part)) where writes(particle_region) do
--TODO: NYI in DSL.
--2D shape of 9 particles, places at 0, 1, 2 in x and y dimension
--zero_core_part(particle_region)
--zero_neighbour_part(particle_region)
--fill(particle_region.interactions, 0)
zero_part_task(particle_region)
c.legion_runtime_issue_execution_fence(__runtime(), __context())

particle_region[0].core_part_space.pos_x = 0.0
particle_region[0].core_part_space.pos_y = 0.0
particle_region[0].core_part_space.cutoff = 1.75

particle_region[1].core_part_space.pos_x = 1.0
particle_region[1].core_part_space.pos_y = 0.0
particle_region[1].core_part_space.cutoff = 1.75

particle_region[2].core_part_space.pos_x = 2.0
particle_region[2].core_part_space.pos_y = 0.0
particle_region[2].core_part_space.cutoff = 1.75


particle_region[3].core_part_space.pos_x = 0.0
particle_region[3].core_part_space.pos_y = 1.0
particle_region[3].core_part_space.cutoff = 1.75

particle_region[4].core_part_space.pos_x = 1.0
particle_region[4].core_part_space.pos_y = 1.0
particle_region[4].core_part_space.cutoff = 1.75

particle_region[5].core_part_space.pos_x = 2.0
particle_region[5].core_part_space.pos_y = 1.0
particle_region[5].core_part_space.cutoff = 1.75


particle_region[6].core_part_space.pos_x = 0.0
particle_region[6].core_part_space.pos_y = 2.0
particle_region[6].core_part_space.cutoff = 1.75

particle_region[7].core_part_space.pos_x = 1.0
particle_region[7].core_part_space.pos_y = 2.0
particle_region[7].core_part_space.cutoff = 1.75

particle_region[8].core_part_space.pos_x = 2.0
particle_region[8].core_part_space.pos_y = 2.0
particle_region[8].core_part_space.cutoff = 1.75



particle_region[0].core_part_space.id = 0
particle_region[1].core_part_space.id = 1
particle_region[2].core_part_space.id = 2
particle_region[3].core_part_space.id = 3
particle_region[4].core_part_space.id = 4
particle_region[5].core_part_space.id = 5
particle_region[6].core_part_space.id = 6
particle_region[7].core_part_space.id = 7
particle_region[8].core_part_space.id = 8
end

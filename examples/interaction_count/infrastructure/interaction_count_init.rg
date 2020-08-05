import "regent"

require("defaults")
require("src/neighbour_search/cell_pair/cell")
require("src/particles/init_part")

local c = regentlib.c

--TODO: This is hopefully what the initialisation function looks like for now
--It can of course read from files instead when necessary.
function initialisation_function(particle_array, space) 

local init_string = rquote
  var particles_space = ispace(int1d, 9)
  var [particle_array] = region(particles_space, part)
  var [space] = region(ispace(int1d, 1), space_config)
  fill([space].{dim_x, dim_y, dim_z}, 0.0)
  
  
  init_space(3.0, 3.0, 3.0, [space])
  particle_initialisation([particle_array])
end
return init_string
end

--TODO: This is hopefully what the finalisation function looks like for now.
--Has loads of excess STUFF still to test this example carefully
function finalisation_function(particle_array, space)

local final_string = rquote
  c.legion_runtime_issue_execution_fence(__runtime(), __context())
  format.println("{}", [particle_array][0].interactions)
  for point in [particle_array].ispace do
    regentlib.assert([particle_array][point].interactions == 8, "test failed")
  end
  
  write_hdf5_snapshot("examples/interaction_count/basic_test.hdf5", [particle_array], [space])
  c.legion_runtime_issue_execution_fence(__runtime(), __context())
  var read_count = read_particle_count("examples/interaction_count/basic_test.hdf5")
  var copy_region = region(ispace(int1d, read_count), part)
  read_hdf5_snapshot("examples/interaction_count/basic_test.hdf5", copy_region)
  
  for point in copy_region do
    regentlib.assert(copy_region[point].core_part_space.pos_x == [particle_array][point].core_part_space.pos_x, "failed on pos_x")
    regentlib.assert(copy_region[point].core_part_space.pos_y == [particle_array][point].core_part_space.pos_y, "failed on pos_y")
    regentlib.assert(copy_region[point].core_part_space.pos_z == [particle_array][point].core_part_space.pos_z, "failed on pos_z")
    regentlib.assert(copy_region[point].core_part_space.cutoff == [particle_array][point].core_part_space.cutoff, "failed on cutoff")
  end
end
return final_string
end


task particle_initialisation(particle_region : region(ispace(int1d), part)) where writes(particle_region) do
--TODO: NYI in DSL.
--2D shape of 9 particles, places at 0, 1, 2 in x and y dimension
zero_core_part(particle_region)
zero_neighbour_part(particle_region)
fill(particle_region.interactions, 0)
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
end

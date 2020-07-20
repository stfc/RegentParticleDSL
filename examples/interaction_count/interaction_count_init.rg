import "regent"

require("defaults")
require("src/neighbour_search/cell_pair/cell")
require("src/particles/init_part")

local c = regentlib.c

task particle_initialisation(particle_region : region(ispace(int1d), part)) where writes(particle_region) do
--TODO: NYI in DSL.
--2D shape of 9 particles, places at 0, 1, 2 in x and y dimension
zero_core_part(particle_region)
zero_neighbour_part(particle_region)
fill(particle_region.interactions, 0)
c.legion_runtime_issue_execution_fence(__runtime(), __context())

particle_region[0].core_part_space.pos_x = 0.0
particle_region[0].core_part_space.pos_y = 0.0
particle_region[0].core_part_space.cutoff = 5.0

particle_region[1].core_part_space.pos_x = 1.0
particle_region[1].core_part_space.pos_y = 0.0
particle_region[1].core_part_space.cutoff = 5.0

particle_region[2].core_part_space.pos_x = 2.0
particle_region[2].core_part_space.pos_y = 0.0
particle_region[2].core_part_space.cutoff = 5.0


particle_region[3].core_part_space.pos_x = 0.0
particle_region[3].core_part_space.pos_y = 1.0
particle_region[3].core_part_space.cutoff = 5.0

particle_region[4].core_part_space.pos_x = 1.0
particle_region[4].core_part_space.pos_y = 1.0
particle_region[4].core_part_space.cutoff = 5.0

particle_region[5].core_part_space.pos_x = 2.0
particle_region[5].core_part_space.pos_y = 1.0
particle_region[5].core_part_space.cutoff = 5.0


particle_region[6].core_part_space.pos_x = 0.0
particle_region[6].core_part_space.pos_y = 2.0
particle_region[6].core_part_space.cutoff = 5.0

particle_region[7].core_part_space.pos_x = 1.0
particle_region[7].core_part_space.pos_y = 2.0
particle_region[7].core_part_space.cutoff = 5.0

particle_region[8].core_part_space.pos_x = 2.0
particle_region[8].core_part_space.pos_y = 2.0
particle_region[8].core_part_space.cutoff = 5.0
end

-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"
require("defaults")

task zero_core_part(particle_region : region(ispace(int1d), part)) where writes(particle_region.core_part_space) do
fill(particle_region.core_part_space.pos_x, 0.0)
fill(particle_region.core_part_space.pos_y, 0.0)
fill(particle_region.core_part_space.pos_z, 0.0)
fill(particle_region.core_part_space.vel_x, 0.0)
fill(particle_region.core_part_space.vel_y, 0.0)
fill(particle_region.core_part_space.vel_z, 0.0)
fill(particle_region.core_part_space.mass, 0.0)
fill(particle_region.core_part_space.cutoff, 0.0)
fill(particle_region.core_part_space.id, int1d(0))
end

-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

require("defaults")
require("src/particles/init_part")

local c = regentlib.c


function initialisation_function(variables, part_count, space_x, space_y, space_z)

local init_string = rquote
  regentlib.assert(part_count > 0, "Particle count must be at least 1")
  regentlib.assert(space_x > 0.0, "Space x dimension must be greater than 0")
  regentlib.assert(space_y > 0.0, "Space y dimension must be greater than 0")
  regentlib.assert(space_z > 0.0, "Space z dimension must be greater than 0")
  var particles_space = ispace(int1d, part_count)
  var [variables.particle_array] = region(particles_space, part)
  var [variables.config] = region(ispace(int1d, 1), config_type)
  fill([variables.config].{space.dim_x, space.dim_y, space.dim_z}, 0.0)

  init_space(space_x, space_y, space_z, [variables.config])
  particle_initialisation([variables.particle_array])
end
return init_string
end

task particle_initialisation(particle_region : region(ispace(int1d), part)) where writes(particle_region) do

 [generate_zero_part_quote(particle_region)];
end

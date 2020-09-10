-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"


require("defaults")


function kick_kernel(part, config)
local kernel = rquote
  part.core_part_space.vel_x = part.core_part_space.vel_x + (part.accel_x * config.space.timestep * 0.5)
  part.core_part_space.vel_y = part.core_part_space.vel_y + (part.accel_y * config.space.timestep * 0.5)
  part.core_part_space.vel_z = part.core_part_space.vel_z + (part.accel_z * config.space.timestep * 0.5)
end
return kernel
end

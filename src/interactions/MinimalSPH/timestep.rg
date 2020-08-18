import "regent"


require("defaults")


function kick_kernel(part, config)
local kernel = rquote
  part.core_part_space.vel_x = part.core_part_space.vel_x + (part.accel_x * config.space.timestep)
  part.core_part_space.vel_y = part.core_part_space.vel_y + (part.accel_y * config.space.timestep)
  part.core_part_space.vel_z = part.core_part_space.vel_z + (part.accel_z * config.space.timestep)
end
return kernel
end

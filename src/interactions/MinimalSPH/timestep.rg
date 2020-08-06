import "regent"


require("defaults")


function kick_kernel(part, space)
local kernel = rquote
  part.core_part_space.vel_x = part.core_part_space.vel_x + (part.accel_x * space.timestep)
  part.core_part_space.vel_y = part.core_part_space.vel_y + (part.accel_y * space.timestep)
  part.core_part_space.vel_z = part.core_part_space.vel_z + (part.accel_z * space.timestep)
end
return kernel
end

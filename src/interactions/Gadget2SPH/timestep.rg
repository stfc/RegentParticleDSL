import "regent"


require("defaults")


function kick_kernel(part, timestep)
local kernel = rquote
  part.v_full_x = part.v_full_x + (part.accel_x * timestep)
  part.v_full_y = part.v_full_y + (part.accel_y * timestep)
  part.v_full_z = part.v_full_z + (part.accel_z * timestep)
end
return kernel
end

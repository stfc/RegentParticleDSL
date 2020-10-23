-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"


require("defaults")
require("src/interactions/MinimalSPH/constants")
require("src/interactions/MinimalSPH/interactions")
local cmath = terralib.includec("math.h")

local fabsf = regentlib.fabs(double)

function kick1_kernel(part, config)
local kernel = rquote
  part.core_part_space.vel_x = part.core_part_space.vel_x + (part.accel_x * config.space.timestep * 0.5)
  part.core_part_space.vel_y = part.core_part_space.vel_y + (part.accel_y * config.space.timestep * 0.5)
  part.core_part_space.vel_z = part.core_part_space.vel_z + (part.accel_z * config.space.timestep * 0.5)
end
return kernel
end

function kick2_kernel(part, config)
local kernel = rquote
  part.core_part_space.vel_x = part.core_part_space.vel_x + (part.accel_x * config.space.timestep * 0.5)
  part.core_part_space.vel_y = part.core_part_space.vel_y + (part.accel_y * config.space.timestep * 0.5)
  part.core_part_space.vel_z = part.core_part_space.vel_z + (part.accel_z * config.space.timestep * 0.5)

  part.pressure = gas_pressure_from_internal_energy(part.rho, part.u)
  part.soundspeed = gas_soundspeed_from_pressure(part.rho, part.pressure)
  part.v_sig = regentlib.fmax(part.v_sig, 2.0 * part.soundspeed)
end
return kernel
end

function drift_kernel(part, config)
local kernel = rquote
  part.core_part_space.pos_x = part.core_part_space.pos_x + (part.core_part_space.vel_x * config.space.timestep )
  part.core_part_space.pos_y = part.core_part_space.pos_y + (part.core_part_space.vel_y * config.space.timestep )
  part.core_part_space.pos_z = part.core_part_space.pos_z + (part.core_part_space.vel_z * config.space.timestep )

  --Drift extra properties
  part.u = part.u + part.u_dt * config.space.timestep
  var h_inv = 1.0 / part.h
  --Predict smoothing length
  var w1 = part.h_dt * h_inv * config.space.timestep
  part.h = part.h * cmath.exp(w1)
  --Predict density
  var w2 = -kernel_dimension * w1
  part.rho = part.rho * cmath.exp(w2)
  part.pressure = gas_pressure_from_internal_energy(part.rho, part.u)
  part.soundspeed = gas_soundspeed_from_pressure(part.rho, part.pressure)
  part.v_sig = regentlib.fmax(part.v_sig, 2.0 * part.soundspeed)
end
return kernel
end

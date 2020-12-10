import "regent"

local soundspeed = global(float, 221.47)
local soundspeed_sq = global(float, soundspeed:get() * soundspeed:get())
local density_reference = global(float, 1000.0)

terra pressure_from_density( density: float) : float
  return soundspeed_sq * (density - density_reference)
end


function timestep(part, config)

local timestep_quote = rquote
  if(part.since_euler >= 50) then
    part.since_euler = 0
  end
  if part.is_wall == WALL then
    part.a_hydro_x = 0.0
    part.a_hydro_y = 0.0
    part.a_hydro_z = 0.0
  end
    part.a_hydro_x = part.a_hydro_x + part.a_const_x
    part.a_hydro_y = part.a_hydro_y + part.a_const_y
    part.a_hydro_z = part.a_hydro_z + part.a_const_z
  --Ignoring shifting for now
  if(part.since_euler == 0) then
    var temp : float = part.rho
    part.rho_t_minus1 = temp
    part.rho = part.rho + part.drho_dt*config.space.timestep
    part.core_part_space.pos_x = part.core_part_space.pos_x + part.core_part_space.vel_x*config.space.timestep + 0.5*part.a_hydro_x*config.space.timestep*config.space.timestep
    part.core_part_space.pos_y = part.core_part_space.pos_y + part.core_part_space.vel_y*config.space.timestep + 0.5*part.a_hydro_y*config.space.timestep*config.space.timestep
    part.core_part_space.pos_z = part.core_part_space.pos_z + part.core_part_space.vel_z*config.space.timestep + 0.5*part.a_hydro_z*config.space.timestep*config.space.timestep

    temp = part.core_part_space.vel_x
    part.core_part_space.vel_x = part.core_part_space.vel_x + part.a_hydro_x*config.space.timestep
    part.v_minus1_x = temp
    
    temp = part.core_part_space.vel_y
    part.core_part_space.vel_y = part.core_part_space.vel_y + part.a_hydro_y*config.space.timestep
    part.v_minus1_y = temp

    temp = part.core_part_space.vel_z
    part.core_part_space.vel_z = part.core_part_space.vel_z + part.a_hydro_z*config.space.timestep
    part.v_minus1_z = temp
    part.pressure = pressure_from_density(part.rho)
  else
    var temp : float = part.rho
    part.rho_t_minus1 = temp
    part.rho = part.rho + part.drho_dt*config.space.timestep
    part.core_part_space.pos_x = part.core_part_space.pos_x + part.core_part_space.vel_x*config.space.timestep + 0.5*part.a_hydro_x*config.space.timestep*config.space.timestep
    part.core_part_space.pos_y = part.core_part_space.pos_y + part.core_part_space.vel_y*config.space.timestep + 0.5*part.a_hydro_y*config.space.timestep*config.space.timestep
    part.core_part_space.pos_z = part.core_part_space.pos_z + part.core_part_space.vel_z*config.space.timestep + 0.5*part.a_hydro_z*config.space.timestep*config.space.timestep

    temp = part.core_part_space.vel_x
    part.core_part_space.vel_x = part.v_minus1_x + 2.0*part.a_hydro_x*config.space.timestep
    part.v_minus1_x = temp


    temp = part.core_part_space.vel_y
    part.core_part_space.vel_y = part.v_minus1_y + 2.0*part.a_hydro_y*config.space.timestep
    part.v_minus1_y = temp

    temp = part.core_part_space.vel_z
    part.core_part_space.vel_z = part.v_minus1_z + 2.0*part.a_hydro_z*config.space.timestep
    part.v_minus1_z = temp

    part.pressure = pressure_from_density(part.rho)
  end
  part.since_euler = part.since_euler+1
end
return timestep_quote
end

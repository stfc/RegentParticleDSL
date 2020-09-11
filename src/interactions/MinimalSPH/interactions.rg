-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

require("defaults")
require("src/neighbour_search/cell_pair/neighbour_search")
sqrt = regentlib.sqrt(float)
cbrt = regentlib.cbrt(float)

--Viscosity parameters
local const_viscosity_beta = 3.0

--Moved to constants.rg
--Cubic spline kernel
--local kernel_gamma = 1.732051
--local kernel_const = 8.0 / 3.0
----Adiabatic index
--local hydro_gamma = 5.0/3.0
--local hydro_gamma_minus_one = hydro_gamma-1.0

__demand(__inline)
task pow_minus_gamma_minus_one(x : float)
  var icbrt = 1.0 / cbrt(x)
  return icbrt * icbrt  
end

--Ideal gas equation
__demand(__inline)
task gas_pressure_from_internal_energy(density : float, u : float)
  return hydro_gamma_minus_one * u * density
end

__demand(__inline)
task gas_soundspeed_from_pressure(density : float, P : float) 
  var density_inv = 1.0 / density
  return sqrt(hydro_gamma * P * density_inv)
end

__demand(__inline)
task gas_entropy_from_internal_energy(density : float, u : float)
  return hydro_gamma_minus_one * u *  pow_minus_gamma_minus_one(density)
end

--3Dimensions, so x^4
terra pow_dimension_plus_one( val : float)
  var x2 = val * val
  return x2 * x2
end

--TODO: Optimise kernel evaluation (Could use C code?)
--Think this should be cubic spline, 3 Dimension form only.
terra eval_kernel_func( ui : float, wi : &float, wi_dx : &float )
  var kernel_gamma_inv = 1.0 / kernel_gamma
  var q3 = 0.0
  var q2 = 0.0
  var q1 = 0.0
  var q0 = 0.0
  var x = ui * kernel_gamma_inv
  var w :float = 0.0
  var w_dx : float = 0.0
  if ui <= 0.5 then
    q3 = 3. 
    q2 = -3.0
    q1 = 0
    q0 = 0.5
  else
    q3 = -1.0
    q2 = 3.0
    q1 = -3.0
    q0 = 1.0
  end

  w = q3 * x + q2
  w_dx = q3

  w_dx = w_dx * x + w
  w = w * x + q1

  w_dx = w_dx * x + w
  w = w * x + q0

  if (w < 0.0) then w = 0.0 end
  if (w_dx > 0.0) then w_dx = 0.0 end
  w = w * kernel_constant * kernel_gamma_inv_dim 
  w_dx = w_dx * kernel_constant * kernel_gamma_inv_dim_plus_one

  @wi = w
  @wi_dx = w_dx
--  @wi = q3 * x + q2
--  @wi_dx = q3
--  --
--  @wi_dx = @wi_dx * x + @wi
--  @wi = @wi * x + q1
--  --
--  @wi_dx = @wi_dx * x + @wi
--  @wi = @wi * x + q0
end

function reset_density(part1, config)
  local kernel = rquote
     part1.wcount = 0.0
     part1.wcount_dh = 0.0
     part1.rho = 0.0
     part1.rho_dh = 0.0
     part1.div_v = 0.0
     part1.rot_v_x = 0.0 
     part1.rot_v_y = 0.0 
     part1.rot_v_z = 0.0 
  end
  return kernel
end

function reset_acceleration(part1, config)
  local kernel = rquote
    part1.accel_x = 0.0
    part1.accel_y = 0.0
    part1.accel_z = 0.0

    part1.u_dt = 0.0
    part1.h_dt = 0.0
    part1.v_sig = 2.0 * part1.soundspeed 
  end
  return kernel

end

function end_force(part1, config)
  local kernel = rquote
    part1.h_dt = part1.h * hydro_dimension_inv
  end
  return kernel
end

function nonsym_density_kernel(part1, part2, r2)
local hydro_dimension = 3.0
local kernel = rquote

  var mj = part2.core_part_space.mass

  var ir = 1.0 / sqrt(r2)
  var r = r2 * ir

  var hi_inv = 1.0 / part1.h
  var ui = r * hi_inv
  var wi : float
  var wi_dx : float
  eval_kernel_func(ui, &wi, &wi_dx)

  part1.rho = (part1.rho + mj * wi)
  part1.rho_dh = part1.rho_dh - (mj * (hydro_dimension * wi + ui * wi_dx))
  --if(part1.core_part_space.id == int1d(136970)) then
  --  format.println("wi = {} wcount = {} ui = {} h = {} p2 = {} x = {} kgi = {} redo = {}", wi, part1.wcount, ui, part1.h, part2.core_part_space.id, ui / kernel_gamma, 1.0 / kernel_gamma, part1.cutoff_update_space.redo)
--  end
  part1.wcount = part1.wcount + wi
--  if(part1.core_part_space.id == int1d(136970)) then
--    format.println("wi = {} wcount = {}", wi, part1.wcount)
--  end
  part1.wcount_dh = part1.wcount_dh - (hydro_dimension * wi + ui * wi_dx)

  var faci = mj * wi_dx * ir

  --Compute dvdr
  var dv_x = part1.core_part_space.vel_x - part2.core_part_space.vel_x
  var dv_y = part1.core_part_space.vel_y - part2.core_part_space.vel_y
  var dv_z = part1.core_part_space.vel_z - part2.core_part_space.vel_z
  var dx_x = part1.core_part_space.pos_x - part2.core_part_space.pos_x
  var dx_y = part1.core_part_space.pos_y - part2.core_part_space.pos_y
  var dx_z = part1.core_part_space.pos_z - part2.core_part_space.pos_z
  var dvdr = dv_x * dx_x + dv_y * dx_y + dv_z * dx_z

  part1.div_v = part1.div_v - (faci * dvdr)

  --Compute dv cross r
  var curlvr_x = dv_y * dx_z - dv_z * dx_y
  var curlvr_y = dv_z * dx_x - dv_x * dx_z
  var curlvr_z = dv_x * dx_y - dv_y * dx_x

  part1.rot_v_x = part1.rot_v_x + (faci * curlvr_x)
  part1.rot_v_y = part1.rot_v_y + (faci * curlvr_y)
  part1.rot_v_z = part1.rot_v_z + (faci * curlvr_z)

end
return kernel
end

function density_kernel(part1, part2, r2)

local hydro_dimension = 3.0
local kernel = rquote
  var mi = part1.core_part_space.mass
  var mj = part2.core_part_space.mass

  var ir = 1.0 / sqrt(r2)
  var r = r2 * ir

  var hi_inv = 1.0 / part1.h
  var ui = r * hi_inv
  var wi : float
  var wi_dx : float
  eval_kernel_func(ui, &wi, &wi_dx)

  part1.rho = (part1.rho + mj * wi)
  part1.rho_dh = part1.rho_dh - (mj * (hydro_dimension * wi + ui * wi_dx))
  part1.wcount = part1.wcount + wi
  part1.wcount_dh = part1.wcount_dh - (hydro_dimension * wi + ui * wi_dx)

  var hj_inv = 1.0 / part2.h
  var uj = r * hj_inv
  var wj : float
  var wj_dx : float
  eval_kernel_func(uj, &wj, &wj_dx)
  part2.rho = part2.rho + mi * wj
  part2.rho_dh = part2.rho_dh - (mi * (hydro_dimension * wj + uj * wj_dx))
  part2.wcount = part2.wcount + wj
  part2.wcount_dh = part2.wcount_dh - (hydro_dimension * wj + uj * wj_dx)

  var faci = mj * wi_dx * ir
  var facj = mi * wj_dx * ir

  --Compute dvdr
  var dv_x = part1.core_part_space.vel_x - part2.core_part_space.vel_x
  var dv_y = part1.core_part_space.vel_y - part2.core_part_space.vel_y
  var dv_z = part1.core_part_space.vel_z - part2.core_part_space.vel_z
  var dx_x = part1.core_part_space.pos_x - part2.core_part_space.pos_x
  var dx_y = part1.core_part_space.pos_y - part2.core_part_space.pos_y
  var dx_z = part1.core_part_space.pos_z - part2.core_part_space.pos_z
  var dvdr = dv_x * dx_x + dv_y * dx_y + dv_z * dx_z

  part1.div_v = part1.div_v - (faci * dvdr)
  part2.div_v = part2.div_v - (facj * dvdr)

  --Compute dv cross r
  var curlvr_x = dv_y * dx_z - dv_z * dx_y
  var curlvr_y = dv_z * dx_x - dv_x * dx_z
  var curlvr_z = dv_x * dx_y - dv_y * dx_x

  part1.rot_v_x = part1.rot_v_x + (faci * curlvr_x)
  part1.rot_v_y = part1.rot_v_y + (faci * curlvr_y)
  part1.rot_v_z = part1.rot_v_z + (faci * curlvr_z)
  
 
  part2.rot_v_x = part2.rot_v_x + (facj * curlvr_x)
  part2.rot_v_y = part2.rot_v_y + (facj * curlvr_y)
  part2.rot_v_z = part2.rot_v_z + (facj * curlvr_z)
end
return kernel
end


function nonsym_force_kernel(part1, part2, r2)

local kernel = rquote

  var mi = part1.core_part_space.mass
  var mj = part2.core_part_space.mass
  var rhoi = part1.rho
  var rhoj = part2.rho
  var pressurei = part1.pressure
  var pressurej = part2.pressure

  var ir = 1.0 / sqrt(r2)
  var r = r2 * ir

--Compute kernel i
  var hi_inv = 1.0 / part1.h
  var hid_inv = pow_dimension_plus_one(hi_inv)
  var ui = r * hi_inv
  var wi : float
  var wi_dx : float
  eval_kernel_func(ui, &wi, &wi_dx)
  var wi_dr = hid_inv * wi_dx

--Compute kernel j
  var hj_inv = 1.0 / part2.h
  var hjd_inv = pow_dimension_plus_one(hj_inv)
  var uj = r * hj_inv
  var wj : float = wj_array[0]
  var wj_dx : float = wj_array[1]
  eval_kernel_func(uj, &wj, &wj_dx)
  var wj_dr = hjd_inv * wj_dx

--h-gradient terms
  var f_i = part1.f
  var f_j = part2.f

--Pressure terms
  var P_over_rho2_i = pressurei / (rhoi * rhoj) * f_i
  var P_over_rho2_j = pressurej / (rhoj * rhoi) * f_j

--Compute soundspeed
  var ci = part1.soundspeed
  var cj = part2.soundspeed

--Compute dvdr
  var dv_x = part1.core_part_space.vel_x - part2.core_part_space.vel_x
  var dv_y = part1.core_part_space.vel_y - part2.core_part_space.vel_y
  var dv_z = part1.core_part_space.vel_z - part2.core_part_space.vel_z
  var dx_x = part1.core_part_space.pos_x - part2.core_part_space.pos_x
  var dx_y = part1.core_part_space.pos_y - part2.core_part_space.pos_y
  var dx_z = part1.core_part_space.pos_z - part2.core_part_space.pos_z
  var dvdr = dv_x * dx_x + dv_y * dx_y + dv_z * dx_z

--Get the balsara term
  var balsara_i = part1.balsara
  var balsara_j = part2.balsara

--Are particles moving towards each other
  var omega_ij = regentlib.fmin(dvdr, 0.0)
  var mu_ij = ir * omega_ij

--Compute the signal velocity
  var v_sig = ci + cj - const_viscosity_beta * mu_ij

--Construct the viscosity term
  var rho_ij = 0.5 * ( rhoi + rhoj)
  var visc = -0.25 * v_sig * mu_ij * (balsara_i + balsara_j) / rho_ij

--Add kernel effect
  var visc_term = 0.5 * visc * (wi_dr + wj_dr) * ir
  var sph_term = (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * ir

  var acceleration = visc_term + sph_term

--Update forces
  part1.accel_x =  part1.accel_x - ( mj * acceleration * dx_x)
  part1.accel_y =  part1.accel_y - ( mj * acceleration * dx_y)
  part1.accel_z =  part1.accel_z - ( mj * acceleration * dx_z)

--Update the signal velocity
  part1.h_dt = part1.h_dt - (mj * dvdr * ir / rhoj * wi_dr)

--Get the time derivative for the internal energy
  var sph_du_term_i = P_over_rho2_i * dvdr * ir * wi_dr
--viscosity term
  var visc_du_term = 0.5 * visc_term * dvdr

--Make the energy equation term
  var du_dt_i = sph_du_term_i + visc_du_term

  part1.u_dt = part1.u_dt + du_dt_i * mj

end
return kernel
end

function force_kernel(part1, part2, r2)

local kernel = rquote

  var mi = part1.core_part_space.mass
  var mj = part2.core_part_space.mass
  var rhoi = part1.rho
  var rhoj = part2.rho
  var pressurei = part1.pressure
  var pressurej = part2.pressure

  var ir = 1.0 / sqrt(r2)
  var r = r2 * ir

--Compute kernel i
  var hi_inv = 1.0 / part1.h
  var hid_inv = pow_dimension_plus_one(hi_inv)
  var ui = r * hi_inv
  var wi : float
  var wi_dx : float
  eval_kernel_func(ui, &wi, &wi_dx)
  var wi_dr = hid_inv * wi_dx

--Compute kernel j
  var hj_inv = 1.0 / part2.h
  var hjd_inv = pow_dimension_plus_one(hj_inv)
  var uj = r * hj_inv
  var wj : float
  var wj_dx : float
  eval_kernel_func(uj, &wj, &wj_dx)
  var wj_dr = hjd_inv * wj_dx

--h-gradient terms
  var f_i = part1.f
  var f_j = part2.f

--Pressure terms
  var P_over_rho2_i = pressurei / (rhoi * rhoj) * f_i
  var P_over_rho2_j = pressurej / (rhoj * rhoi) * f_j

--Compute soundspeed
  var ci = part1.soundspeed
  var cj = part2.soundspeed

--Compute dvdr
  var dv_x = part1.core_part_space.vel_x - part2.core_part_space.vel_x
  var dv_y = part1.core_part_space.vel_y - part2.core_part_space.vel_y
  var dv_z = part1.core_part_space.vel_z - part2.core_part_space.vel_z
  var dx_x = part1.core_part_space.pos_x - part2.core_part_space.pos_x
  var dx_y = part1.core_part_space.pos_y - part2.core_part_space.pos_y
  var dx_z = part1.core_part_space.pos_z - part2.core_part_space.pos_z
  var dvdr = dv_x * dx_x + dv_y * dx_y + dv_z * dx_z

--Get the balsara term
  var balsara_i = part1.balsara
  var balsara_j = part2.balsara

--Are particles moving towards each other
  var omega_ij = regentlib.fmin(dvdr, 0.0)
  var mu_ij = ir * omega_ij

--Compute the signal velocity
  var v_sig = ci + cj - const_viscosity_beta * mu_ij

--Construct the viscosity term
  var rho_ij = 0.5 * ( rhoi + rhoj)
  var visc = -0.25 * v_sig * mu_ij * (balsara_i + balsara_j) / rho_ij
  
--Add kernel effect
  var visc_term = 0.5 * visc * (wi_dr + wj_dr) * ir
  var sph_term = (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * ir

  var acceleration = visc_term + sph_term

--Update forces
  part1.accel_x =  part1.accel_x - ( mj * acceleration * dx_x)
  part1.accel_y =  part1.accel_y - ( mj * acceleration * dx_y)
  part1.accel_z =  part1.accel_z - ( mj * acceleration * dx_z)
  
  part2.accel_x =  part2.accel_x + ( mi * acceleration * dx_x)
  part2.accel_y =  part2.accel_y + ( mi * acceleration * dx_y)
  part2.accel_z =  part2.accel_z + ( mi * acceleration * dx_z)

--Update the signal velocity
  part1.h_dt = part1.h_dt - (mj * dvdr * ir / rhoj * wi_dr)
  part2.h_dt = part2.h_dt - (mi * dvdr * ir / rhoi * wj_dr)

--Get the time derivative for the internal energy
  var sph_du_term_i = P_over_rho2_i * dvdr * ir * wi_dr
  var sph_du_term_j = P_over_rho2_j * dvdr * ir * wj_dr
--viscosity term
  var visc_du_term = 0.5 * visc_term * dvdr

--Make the energy equation term
  var du_dt_i = sph_du_term_i + visc_du_term
  var du_dt_j = sph_du_term_j+ visc_du_term

  part1.u_dt = part1.u_dt + du_dt_i * mj
  part2.u_dt = part2.u_dt + du_dt_j * mi
 
 
end
return kernel
end

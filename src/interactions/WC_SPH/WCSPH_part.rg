-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"
require("src/particles/core_part")



fspace part{
  core_part_space : core_part,
  neighbour_part_space : neighbour_part,
  --smoothing length
  h : float,
  --internal energy
  u : float,
  --Density
  rho : float,
  rho_base : float,
  rho_t_minus1 : float,
  viscosity : float, 
  pressure : float,
  soundspeed : float,
  --Stress tensor
  tau_xx : float,
  tau_xy : float,
  tau_xz : float,
  tau_yy : float,
  tau_yz : float,
  tau_zz : float,
  --Velocity gradient
  grad_v_xx : float,
  grad_v_xy : float,
  grad_v_xz : float,
  grad_v_yy : float,
  grad_v_yz : float,
  grad_v_zz : float,
  --Particle type
  is_wall : int32,
  since_euler : int32,
  max_visc : float,
  --Acceleration
  a_hydro_x : float,
  a_hydro_y : float,
  a_hydro_z : float,
  --Previous timestep data
  v_minus1_x : float,
  v_minus1_y : float,
  v_minus1_z : float,

}

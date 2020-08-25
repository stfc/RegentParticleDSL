import "regent"

require("src/particles/core_part")
require("src/neighbour_search/cell_pair/import_cell_pair")

--Extra data required for the cutoff update
fspace cutoff_update{
  redo : int1d,
  left : float,
  right : float,
  h_0 : float
}


fspace part{
  neighbour_part_space : neighbour_part,
--Contains position, velocity and mass
  core_part_space : core_part,
  cutoff_update_space : cutoff_update,
  --SPH h - used to compute the cutoff radius since the cutoff
  --in the core part space is technically the search radius
  h : float,
  --Particle acceleration
  accel_x : float,
  accel_y : float,
  accel_z : float,
  --internal energy
  u : float,
  u_dt : float,
  --density
  rho : float,
  --weighted neighbour count and derivative
  wcount : float,
  wcount_dh : float,
  --derivate of density w.r.t to h
  rho_dh : float,
  --Velocity divergence
  div_v : float,
  --velocity curl
  rot_v_x : float,
  rot_v_y : float,
  rot_v_z : float,
  --Balsara switch
  balsara : float,
  --"Grad h term"
  f : float,
  --Pressure
  pressure : float,
  --Particle soundspeed
  soundspeed : float,
  --Signal velocity
  v_sig : float,
  --Time derivative of the smoothing length
  h_dt : float,
}

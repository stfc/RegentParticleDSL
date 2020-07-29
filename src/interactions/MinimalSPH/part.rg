import "regent"

require("src/particles/core_part")
require("src/neighbour_search/cell_pair/import_cell_pair")

--Extra data required for the cutoff update
fspace cutoff_update{
  redo : int1d,
  left : double,
  right : double,
  h_0 : double
}


fspace part{
  neighbour_part_space : neighbour_part,
--Contains position, velocity and mass
  core_part_space : core_part,
  cutoff_update_space : cutoff_update,
  --SPH h - used to compute the cutoff radius since the cutoff
  --in the core part space is technically the search radius
  h : double,
  --Particle acceleration
  accel_x : double,
  accel_y : double,
  accel_z : double,
  --internal energy
  u : double,
  u_dt : double,
  --density
  rho : double,
  --weighted neighbour count and derivative
  wcount : double,
  wcount_dh : double,
  --derivate of density w.r.t to h
  rho_dh : double,
  --Velocity divergence
  div_v : double,
  --velocity curl
  rot_v_x : double,
  rot_v_y : double,
  rot_v_z : double,
  --Balsara switch
  balsara : double,
  --"Grad h term"
  f : double,
  --Pressure
  pressure : double,
  --Particle soundspeed
  soundspeed : double,
  --Signal velocity
  v_sig : double,
  --Time derivative of the smoothing length
  h_dt : double,
}

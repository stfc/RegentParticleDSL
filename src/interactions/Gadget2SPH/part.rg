import "regent"

require("src/particles/core_part")
require("src/neighbour_search/cell_pair/import_cell_pair")

fspace part{
  neighbour_part_space : neighbour_part,
  core_part_space : core_part,
  --velocity at the last full step
  v_full_x : double,
  v_full_y : double,
  v_full_z : double,
  --Particle position
  pos_x : double,
  pos_y : double,
  pos_z : double,
  --Predicted velocity
  vel_x : double,
  vel_y : double,
  vel_z : double,
  --Particle acceleration
  accel_x : double,
  accel_y : double,
  accel_z : double,
  --Particle mass
  mass : double,
  --density
  rho : double,
  --weighted neighbour count and derivative
  wcount : double,
  wcount_dh : double,
  --derivate of density w.r.t to h
  rho_dh : double,
  --velocity curl
  rot_v_x : double,
  rot_v_y : double,
  rot_v_z : double,
  --velocity divergence
  div_v : double,
  --Balsara switch
  balsara : double,
  --"Grad h term"
  f : double,
  --Pressure over density squared
  P_over_rho2 : double,
  --Particle soundspeed
  soundspeed : double,
  --Signal velocity
  v_sig : double,
  --Time derivative of the smoothing length
  h_dt : double
}

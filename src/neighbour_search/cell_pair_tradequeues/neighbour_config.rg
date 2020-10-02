-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

fspace neighbour_config_type{
  cell_dim_x : double,
  cell_dim_y : double,
  cell_dim_z : double,
  x_cells : int,
  y_cells : int, 
  z_cells : int,
  max_cutoff : double
}

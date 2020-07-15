import "regent"

--This file contains the core particle "type", which all particles used
--must expand from.


--All particles must have mass, position, velocity and ID
fspace core_part{
  pos_x : double,
  pos_y : double,
  pos_z : double,
  mass : double,
  vel_x : double,
  vel_y : double,
  vel_z : double,
  cutoff : double,
  id : int1d
}

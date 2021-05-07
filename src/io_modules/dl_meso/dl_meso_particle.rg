-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"
require("src/particles/core_part")


fspace part{
  core_part_space : core_part,
  neighbour_part_space : neighbour_part,

  lab : int32, -- This could be jsut core_part_space.id - we'll see
  ltp : int, --Particle species
  ltm : int, --Particle molecule type (TODO: NYI)
  atmnam : int8[25],
 
  fxx : double, --Force components
  fyy : double,
  fzz : double,

  fvx : double, --variable particle force
  fvy : double,
  fvz : double,
}

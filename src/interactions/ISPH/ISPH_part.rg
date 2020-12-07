-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"
require("src/particles/core_part")

print("ISPH Part imported")

fspace part{
  core_part_space : core_part,
  neighbour_part_space : neighbour_part,
  pressure : double,
  volume : double,
  divergence : double,
  interactions : int
}

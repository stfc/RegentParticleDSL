-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

require("src/particles/core_part")
require("src/neighbour_search/cell_pair/import_cell_pair")

fspace part{
  neighbour_part_space : neighbour_part,
  core_part_space : core_part,
  --Only extra variable we need is a counter of how many interactions we have!
  interactions : uint32
}

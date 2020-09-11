-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

--require("src/particles/core_part")
--require("src/neighbour_search/default_neighbour_search")
require("defaults")

--Definition of the default particle type with no extra data.
fspace part{
--We are required to include the neighbour_part and core_part types in our particle
--declaration
  neighbour_part_space : neighbour_part,
  core_part_space : core_part,
}

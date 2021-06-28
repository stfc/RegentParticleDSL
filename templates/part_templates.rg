-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

--This file shows an example construction of a custom particle type, used in a computation


--This structure must always be called "part"
fspace part{
--We are required to include the neighbour_part and core_part types in our particle
--declaration
  neighbour_part_space : neighbour_part,
  core_part_space : core_part,
--Any extra defitions we want for our computation can also be added, e.g.
  extra_variable_1 : double,
  extra_variable_2 : float,
  extra_variable_3 : int32

}

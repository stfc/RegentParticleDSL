-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

--This file contains the neighbour_part type definition for the
--cell pair neighbour method. All particles used for computation
--must include neighbour_part

--We need to store information saying which cell each particle is
--contained in
fspace neighbour_part{

 cell_id : int3d

}

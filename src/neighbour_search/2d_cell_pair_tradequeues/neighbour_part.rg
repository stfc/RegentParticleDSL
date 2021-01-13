-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

--This file contains the neighbour_part type definition for the
--cell pair neighbour method. All particles used for computation
--must include neighbour_part
--print("Included the 2D part")

--We need to store information saying which cell each particle is
--contained in
fspace neighbour_part{

 cell_id : int2d,
 _valid : bool,
--FIXME: Think about adding the following fields to enable optimisations
-- old_pos_x : double,
-- old_pos_y : double,
-- old_pos_z : double,

}

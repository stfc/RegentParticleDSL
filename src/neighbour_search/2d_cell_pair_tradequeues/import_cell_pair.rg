-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

--To use cell_pair neighbour finding you must include this file in
--the require decelarations.
require("src/neighbour_search/2d_cell_pair_tradequeues/neighbour_part")
require("src/neighbour_search/2d_cell_pair_tradequeues/neighbour_config")
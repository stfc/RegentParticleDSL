-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

DSL_DIMENSIONALITY = 2
--To use cell_pair neighbour finding you must include this file in
--the require decelarations.
require("src/neighbour_search/cell_pair_tradequeues_nonperiod/2d_neighbour_part")
require("src/neighbour_search/cell_pair_tradequeues_nonperiod/2d_neighbour_config")

neighbour_search_validity = terralib.newlist()
neighbour_search_validity:insert({field="neighbour_part_space._valid", result=true})

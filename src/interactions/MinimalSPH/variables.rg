-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

require("defaults")

local variables = {}

variables.config = regentlib.newsymbol("config")
variables.particle_array = regentlib.newsymbol("particle_region")
variables.cell_partition = regentlib.newsymbol("cell_partition")
variables.cell_space = regentlib.newsymbol("cell_space")
variables.cell_partition2 = regentlib.newsymbol("cell_partition2")

return variables

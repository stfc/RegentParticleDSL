import "regent"

require("defaults")
require("src/interactions/MinimalSPH/interactions")
require("src/interactions/MinimalSPH/tasks")
require("src/interactions/MinimalSPH/hdf5_io")
require("src/interactions/MinimalSPH/timestep")

local variables = require("src/interactions/MinimalSPH/variables")
local format = require("std/format")

task main()
[initialisation("/home/aidan/swiftsim/examples/HydroTests/SodShock_3D/sodShock.hdf5", variables.particle_array, variables.space)]
format.println("{}", [variables.particle_array][0].core_part_space.id)
end

regentlib.start(main)

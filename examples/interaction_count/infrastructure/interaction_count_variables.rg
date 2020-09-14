import "regent"

local variables = {}

variables.config = regentlib.newsymbol("config")
variables.particle_array = regentlib.newsymbol("particle_region")
variables.cell_space = regentlib.newsymbol("cell_space")
return variables

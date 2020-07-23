import "regent"

require("defaults")
--require("src/neighbour_search/cell_pair/import_cell_pair")
require("src/neighbour_search/cell_pair/neighbour_search")
require("src/neighbour_search/cell_pair/cell")
require("src/interactions/Gadget2SPH/interactions")
require("src/interactions/Gadget2SPH/timestep")

local density_symmetric_task = generate_symmetric_pairwise_task(iact_density)

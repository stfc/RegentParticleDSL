import "regent"

require("defaults")
require("src/neighbour_search/cell_pair/import_cell_pair")
require("src/neighbour_search/cell_pair/neighbour_search")
require("examples/interaction_count/interaction_count_kernel")

task particle_initialisation()
--TODO: NYI in DSL or here.
end

local count_interaction_task = generate_symmetric_pairwise_task( symmetric_interaction_count_kernel )
local interaction_tasks_runner = run_symmetric_pairwise_task( count_interaction_task )

task main_task()

--We will use 9 particles, with a large cutoff. Various of these functions will not be user-implemented later.
var particles_space = ispace(int1d, 9)
var particle_region = region(particles_space, part)

end

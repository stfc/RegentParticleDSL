import "regent"

require("defaults")
require("src/neighbour_search/cell_pair/import_cell_pair")
require("src/neighbour_search/cell_pair/neighbour_search")
require("src/neighbour_search/cell_pair/cell")
require("examples/interaction_count/interaction_count_kernel")
require("examples/interaction_count/interaction_count_init")

local c = regentlib.c
format = require("std/format")

--Create the tasks and runner from the kernel. This will become a single function call soon.
local count_interaction_task = generate_symmetric_pairwise_task( symmetric_interaction_count_kernel )
local interaction_tasks_runner = run_symmetric_pairwise_task( count_interaction_task )

task main_task()

--We will use 9 particles, with a large cutoff. Various of these functions will not be user-implemented later.
var particles_space = ispace(int1d, 9)
var particle_region = region(particles_space, part)

particle_initialisation(particle_region)

--We will set cell size to be 1.0 for now. Not periodic or anything for now so no idea of global cell size for now.
--NYI: The task/library should choose the cell size itself.
particles_to_cell_launcher(particle_region, 1.0, 1.0, 1.0)
--Generate the cell partition. This ideally will not be done by the user, or will be "hidden"
var cell_partition = update_cell_partitions(particle_region, 3, 3, 3)

--Run the interaction tasks "timestep". We could do this multiple types and know it will be correct due to 
--the programming model
interaction_tasks_runner(particle_region, cell_partition)

c.legion_runtime_issue_execution_fence(__runtime(), __context())
format.println("{}", particle_region[0].interactions)
regentlib.assert(particle_region[0].interactions == 8, "test failed")
regentlib.assert(particle_region[6].interactions == 8, "test failed")

end

regentlib.start(main_task)

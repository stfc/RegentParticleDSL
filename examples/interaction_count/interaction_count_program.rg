import "regent"

require("defaults")
require("src/neighbour_search/cell_pair/import_cell_pair")
require("src/neighbour_search/cell_pair/neighbour_search")
require("src/neighbour_search/cell_pair/cell")
require("examples/interaction_count/interaction_count_kernel")
require("examples/interaction_count/infrastructure/interaction_count_init")
require("examples/interaction_count/infrastructure/interaction_count_IO")

local c = regentlib.c
format = require("std/format")
variables = require("examples/interaction_count/infrastructure/interaction_count_variables")
local stdlib = terralib.includec("stdlib.h")

--Create the tasks and runner from the kernel. You can choose to use symmetric or asymmetric as desired here.
--local interaction_tasks_runner = create_symmetric_pairwise_runner( symmetric_interaction_count_kernel )
local interaction_tasks_runner = create_asymmetric_pairwise_runner( asymmetric_interaction_count_kernel )

task main_task()

  --We will use 9 particles, with a large cutoff. The initialisation is currently declared in the interaction_count_init file.
  --In the future much of this will be abstracted into the DSL, though user-specified particles can be done easily enough, and
  --example code to do this will be added.
  [initialisation_function(variables.particle_array, variables.config)];
  
  --We will set cell size to be 1.0 for now. Not periodic or anything for now so no idea of global cell size for now.
  --TODO: NYI: The DSL library should choose the cell size itself.
  --This code will be abstracted into a new function on a per-neighbour finding system later, but with the same structure.
  particles_to_cell_launcher(variables.particle_array, variables.config);
  --Generate the cell partition. This ideally will not be done by the user, or will be "hidden"
  var cell_partition = update_cell_partitions(variables.particle_array, variables.config);
  
  --Run the interaction tasks "timestep". We could do this multiple types and know it will be correct due to 
  --the programming model
  interaction_tasks_runner(variables.particle_array, cell_partition, variables.config);

  --This finalisation function is declared in the interaction_count_init file  
  [finalisation_function(variables.particle_array, variables.config)];
end

regentlib.start(main_task)

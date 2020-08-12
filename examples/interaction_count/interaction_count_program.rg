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

--Create the tasks and runner from the kernel. You can choose to use symmetric or asymmetric as desired here.
--local interaction_tasks_runner = create_symmetric_pairtask_runner( symmetric_interaction_count_kernel )
local interaction_tasks_runner = create_asymmetric_pairwise_runner( asymmetric_interaction_count_kernel )

task main_task()

  --We will use 9 particles, with a large cutoff. Various of these functions will not be user-implemented later.
  --TODO: These lines should become a single call to initialisation, however having issues with the metaprogramming
  --to get this to work right now
  var particles_space = ispace(int1d, 9)
  var [variables.particle_array] = region(particles_space, part)
  var [variables.config] = region(ispace(int1d, 1), config_type)
  fill([variables.config].{space.dim_x, space.dim_y, space.dim_z}, 0.0)
  init_space(3.0, 3.0, 3.0, [variables.config])
  particle_initialisation([variables.particle_array])
  --We will set cell size to be 1.5 for now. 
  --TODO: NYI: The task/library should choose the cell size itself.
 -- var cell_count = compute_cell_count([variables.config], 1.5, 1.5, 1.5)
initialise_cells(variables.config , variables.particle_array)
  particles_to_cell_launcher( variables.particle_array, variables.config)
  --[initialisation_function(variables.particle_array, variables.space)]
  --[initialise]

  ---------------------------------
  --END OF BOILERPLATE INIT CODE --
  ---------------------------------
  
  --Generate the cell partition. This ideally will not be done by the user, or will be "hidden"
  var cell_partition = update_cell_partitions([variables.particle_array], [variables.config])
  
  --Run the interaction tasks "timestep". We could do this multiple types and know it will be correct due to 
  --the programming model
  interaction_tasks_runner([variables.particle_array], cell_partition, [variables.config])
  
  -----------------------------------
  --START OF BOIILERPLATE FIN CODE --
  -----------------------------------
  --TODO: This will be replaced with a finalisation call, however having issues with the metaprogramming
  --to get this to work right now
  --[finalisation(variables.particle_array, variables.space)]
  c.legion_runtime_issue_execution_fence(__runtime(), __context())
  format.println("{}", [variables.particle_array][0].interactions)
  for point in [variables.particle_array].ispace do
    regentlib.assert([variables.particle_array][point].interactions == 8, "test failed")
  end
  
  write_hdf5_snapshot("examples/interaction_count/basic_test.hdf5", [variables.particle_array], [variables.config])
  c.legion_runtime_issue_execution_fence(__runtime(), __context())
  var read_count = read_particle_count("examples/interaction_count/basic_test.hdf5")
  var copy_region = region(ispace(int1d, read_count), part)
  read_hdf5_snapshot("examples/interaction_count/basic_test.hdf5", copy_region)
  
  for point in copy_region do
  regentlib.assert(copy_region[point].core_part_space.pos_x == [variables.particle_array][point].core_part_space.pos_x, "failed on pos_x")
  regentlib.assert(copy_region[point].core_part_space.pos_y == [variables.particle_array][point].core_part_space.pos_y, "failed on pos_y")
  regentlib.assert(copy_region[point].core_part_space.pos_z == [variables.particle_array][point].core_part_space.pos_z, "failed on pos_z")
  regentlib.assert(copy_region[point].core_part_space.cutoff == [variables.particle_array][point].core_part_space.cutoff, "failed on cutoff")
  end
end

regentlib.start(main_task)

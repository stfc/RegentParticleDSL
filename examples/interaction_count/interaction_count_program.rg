import "regent"

require("defaults")
require("src/neighbour_search/cell_pair/import_cell_pair")
require("src/neighbour_search/cell_pair/neighbour_search")
require("src/neighbour_search/cell_pair/cell")
require("examples/interaction_count/interaction_count_kernel")
require("examples/interaction_count/interaction_count_init")
require("examples/interaction_count/interaction_count_IO")

local c = regentlib.c
format = require("std/format")

--Create the tasks and runner from the kernel. This will become a single function call soon.
--local count_interaction_task = generate_symmetric_pairwise_task( symmetric_interaction_count_kernel )
--local interaction_tasks_runner = run_symmetric_pairwise_task( count_interaction_task )
--local interaction_tasks_runner = create_symmetric_pairtask_runner( symmetric_interaction_count_kernel )
local interaction_tasks_runner = create_asymmetric_pairwise_runner( asymmetric_interaction_count_kernel )

task main_task()

--We will use 9 particles, with a large cutoff. Various of these functions will not be user-implemented later.
var particles_space = ispace(int1d, 9)
var particle_region = region(particles_space, part)
var space = region(ispace(int1d, 1), space_config)
fill(space.{dim_x, dim_y, dim_z}, 0.0)


init_space(3.0, 3.0, 3.0, space)
particle_initialisation(particle_region)

--We will set cell size to be 1.0 for now. Not periodic or anything for now so no idea of global cell size for now.
--NYI: The task/library should choose the cell size itself.
particles_to_cell_launcher(particle_region, 1.5, 1.5, 1.5)
--Generate the cell partition. This ideally will not be done by the user, or will be "hidden"
var cell_partition = update_cell_partitions(particle_region, 2, 2, 2)

--Run the interaction tasks "timestep". We could do this multiple types and know it will be correct due to 
--the programming model
interaction_tasks_runner(particle_region, cell_partition, space)

c.legion_runtime_issue_execution_fence(__runtime(), __context())
format.println("{}", particle_region[0].interactions)
for point in particle_region.ispace do
  regentlib.assert(particle_region[point].interactions == 8, "test failed")
end
--regentlib.assert(particle_region[6].interactions == 8, "test failed")

write_hdf5_snapshot("examples/interaction_count/basic_test.hdf5", particle_region, space)
c.legion_runtime_issue_execution_fence(__runtime(), __context())
--var copy_region = read_hdf5_snapshot("test.hdf5")
var read_count = read_particle_count("examples/interaction_count/basic_test.hdf5")
var copy_region = region(ispace(int1d, read_count), part)
read_hdf5_snapshot("examples/interaction_count/basic_test.hdf5", copy_region)

for point in copy_region do
regentlib.assert(copy_region[point].core_part_space.pos_x == particle_region[point].core_part_space.pos_x, "failed on pos_x")
regentlib.assert(copy_region[point].core_part_space.pos_y == particle_region[point].core_part_space.pos_y, "failed on pos_y")
regentlib.assert(copy_region[point].core_part_space.pos_z == particle_region[point].core_part_space.pos_z, "failed on pos_z")
regentlib.assert(copy_region[point].core_part_space.cutoff == particle_region[point].core_part_space.cutoff, "failed on cutoff")
end
end

regentlib.start(main_task)

import "regent"

require("defaults")
require("src/interactions/MinimalSPH/interactions")
require("src/interactions/MinimalSPH/tasks")
require("src/interactions/MinimalSPH/hdf5_io")
require("src/interactions/MinimalSPH/timestep")

local variables = require("src/interactions/MinimalSPH/variables")
local format = require("std/format")
local c = regentlib.c

local density_task = create_asymmetric_pairwise_runner(nonsym_density_kernel)
local timestep_task = run_per_particle_task( kick_kernel )

__forbid(__inline)
task say_hello(time : float)
  format.println("HELLO TIME {}", time)
end

task main()
--[initialisation("/home/aidan/swiftsim/examples/HydroTests/SodShock_3D/sodShock.hdf5", variables.particle_array, variables.space)]
  var filename = "/home/aidan/swiftsim/examples/HydroTests/SodShock_3D/sodShock.hdf5"
  var count = read_particle_count(filename)
  format.println("Initialising SPH from {} with {} hydro particles", filename, count)
  var particles_space = ispace(int1d, count)
  var [variables.particle_array] = region(particles_space, part)
  var [variables.space] = region(ispace(int1d, 1), space_config)
  fill([variables.space].{dim_x, dim_y, dim_z, timestep}, 0.0)
  read_hdf5_snapshot(filename, count, [variables.particle_array], [variables.space])

  format.println("{} {} {}", [variables.space][0].dim_x, [variables.space][0].dim_y, [variables.space][0].dim_z)
  --Make 5x5x5 cells for now (chosen arbitrarily). NB This wouldn't not be sensible for optimised version, but the neighbour search will
  --be involved in cell size choices for real cases
  particles_to_cell_launcher([variables.particle_array],  [variables.space][0].dim_x/5.0,  [variables.space][0].dim_y/5.0,  [variables.space][0].dim_z/5.0)

  var time : double = 0.0
  var endtime : double = 1.0
  [variables.space][0].timestep = 0.01
  while time < endtime do
    var cell_partition = update_cell_partitions([variables.particle_array], 5, 5, 5)
    --density_task([variables.particle_array], cell_partition, [variables.space])
    timestep_task([variables.particle_array], cell_partition, [variables.space])
    c.legion_runtime_issue_execution_fence(__runtime(), __context())
    say_hello(time)
    time = time + 0.01
--    __delete(cell_partition)
  end  
  
end

regentlib.start(main)

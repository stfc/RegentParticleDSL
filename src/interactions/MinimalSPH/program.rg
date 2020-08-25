import "regent"

require("defaults")
require("src/interactions/MinimalSPH/interactions")
require("src/interactions/MinimalSPH/tasks")
require("src/interactions/MinimalSPH/hdf5_io")
require("src/interactions/MinimalSPH/timestep")

local variables = require("src/interactions/MinimalSPH/variables")
local format = require("std/format")
local c = regentlib.c

local density_task = create_asymmetric_pairwise_runner( nonsym_density_kernel )
local force_task = create_symmetric_pairwise_runner( force_kernel )
local timestep_task = run_per_particle_task( kick_kernel )
local reset_density_task = run_per_particle_task( reset_density )
local reset_force_task = run_per_particle_task( reset_acceleration )

__forbid(__inline)
task say_hello(time : float)
  format.println("HELLO TIME {}", time)
end

task get_time() : double

  return c.legion_get_current_time_in_micros()
end

task main()
--[initialisation("/home/aidan/swiftsim/examples/HydroTests/SodShock_3D/sodShock.hdf5", variables.particle_array, variables.space)]
  var filename = "/home/aidan/swiftsim/examples/HydroTests/SodShock_3D/sodShock.hdf5"
  var count = read_particle_count(filename)
  format.println("Initialising SPH from {} with {} hydro particles", filename, count)
  var particles_space = ispace(int1d, count)
  var [variables.particle_array] = region(particles_space, part)
  var [variables.config] = region(ispace(int1d, 1), config_type)
  fill([variables.config].{space.dim_x, space.dim_y, space.dim_z, space.timestep}, 0.0)
  read_hdf5_snapshot(filename, count, [variables.particle_array], [variables.config])

  format.println("{} {} {}", [variables.config][0].space.dim_x, [variables.config][0].space.dim_y, [variables.config][0].space.dim_z)
  --Make 5x5x5 cells for now (chosen arbitrarily). NB This wouldn't not be sensible for optimised version, but the neighbour search will
  --be involved in cell size choices for real cases
initialise_cells(variables.config , variables.particle_array)
  particles_to_cell_launcher( variables.particle_array, variables.config)

  var cell_partition1 = update_cell_partitions([variables.particle_array], [variables.config])
  --Initialisation
  reset_density_task(variables.particle_array, cell_partition1, variables.config)
  reset_force_task(variables.particle_array, cell_partition1, variables.config)

  --Do the zero timestep to setup the IC
  density_task(variables.particle_array, cell_partition1, variables.config)
  update_cutoffs_launcher(variables.particle_array, cell_partition, variables.config)
  var timestep = compute_timestep_launcher(variables.particle_array, cell_partition, variables.config)

  var time : double = 0.0
  var endtime : double = 0.0001
  [variables.config][0].space.timestep = 0.01
  c.legion_runtime_issue_execution_fence(__runtime(), __context())
  var start_time = get_time()
  c.legion_runtime_issue_execution_fence(__runtime(), __context())
  __delete(cell_partition1)
  while time < endtime do
    var cell_partition = update_cell_partitions(variables.particle_array, variables.config)
    reset_density_task(variables.particle_array, cell_partition, variables.config)
    density_task(variables.particle_array, cell_partition, variables.config)
    update_cutoffs_launcher(variables.particle_array, cell_partition, variables.config)
    reset_force_task(variables.particle_array, cell_partition, variables.config)
    force_task( variables.particle_array, cell_partition, variables.config)
  --  timestep_task([variables.particle_array], cell_partition, [variables.space])
    c.legion_runtime_issue_execution_fence(__runtime(), __context())
    say_hello(time)
    time = time + config[0].space.timestep
    config[0].space.timestep = compute_timestep_launcher(variables.particle_array, cell_partition, variables.config)
    if(endtime - time > config[0].space.timestep) then
      config[0].space.timestep = endtime - time
    end
    format.println("timestep is {}", timestep)
    __delete(cell_partition)
  end  
    c.legion_runtime_issue_execution_fence(__runtime(), __context())
  var end_time = get_time()
 
  format.println("Computation took {} seconds.", (end_time - start_time)/1000000.0)
 
end

regentlib.start(main)

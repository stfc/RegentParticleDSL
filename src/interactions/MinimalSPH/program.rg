-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

require("defaults")
require("src/interactions/MinimalSPH/interactions")
require("src/interactions/MinimalSPH/tasks")
require("src/interactions/MinimalSPH/hdf5_io")
require("src/interactions/MinimalSPH/timestep")

local variables = require("src/interactions/MinimalSPH/variables")
local format = require("std/format")
local c = regentlib.c

local density_task = create_asymmetric_pairwise_runner( nonsym_density_kernel, variables.config, variables.cell_partition )
local force_task = create_symmetric_pairwise_runner( force_kernel, variables.config, variables.cell_partition )
local timestep_task = run_per_particle_task( kick_kernel, variables.config, variables.cell_partition )
local kick1_task = run_per_particle_task( kick_kernel, variables.config, variables.cell_partition2 )
local reset_density_task = run_per_particle_task( reset_density, variables.config, variables.cell_partition )
local reset_force_task = run_per_particle_task( reset_acceleration, variables.config, variables.cell_partition )

local initial_density_reset = run_per_particle_task( reset_density, variables.config, variables.cell_space )
local initial_force_reset = run_per_particle_task( reset_density, variables.config, variables.cell_space )
local initial_density_task = create_asymmetric_pairwise_runner( nonsym_density_kernel, variables.config, variables.cell_space )
local initial_timestep_task = run_per_particle_task( kick_kernel, variables.config, variables.cell_space )

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

  var [variables.cell_space] = update_cell_partitions([variables.particle_array], [variables.config])
  --Initialisation
  [initial_density_reset]
  [initial_force_reset]

  --Do the zero timestep to setup the IC
  variables.config[0].space.timestep = 0.00
  [initial_density_task]
  update_cutoffs_launcher(variables.particle_array, variables.cell_space, variables.config)
  [initial_timestep_task]
  variables.config[0].space.timestep = compute_timestep_launcher(variables.particle_array, cell_partition1, variables.config)
  
  var time : double = 0.0
  var endtime : double = 0.02
  c.legion_runtime_issue_execution_fence(__runtime(), __context())
  var start_time = get_time()
  format.println("timestep computed is {}", variables.config[0].space.timestep)
  c.legion_runtime_issue_execution_fence(__runtime(), __context())
  __delete(cell_partition1)
  while time < endtime do
    --TODO Fix update_cell_partitions
    var [variables.cell_partition2] = update_cell_partitions(variables.particle_array, variables.config)
    --first kick
    [kick1_task]
    __delete([variables.cell_partition2])
    var [variables.cell_partition] = update_cell_partitions(variables.particle_array, variables.config)
    [reset_density_task]
    [density_task]
    update_cutoffs_launcher(variables.particle_array, cell_partition, variables.config)
    [reset_force_task]
    [force_task]
    --2nd kick
    [timestep_task]
  --  timestep_task([variables.particle_array], cell_partition, [variables.space])
    c.legion_runtime_issue_execution_fence(__runtime(), __context())
    say_hello(time)
    time = time + variables.config[0].space.timestep
    variables.config[0].space.timestep = compute_timestep_launcher(variables.particle_array, cell_partition, variables.config)
    if(endtime - time > variables.config[0].space.timestep) then
      variables.config[0].space.timestep = endtime - time
    end
    format.println("timestep is {}", variables.config[0].space.timestep)
    __delete([variables.cell_partition])
  end  
    c.legion_runtime_issue_execution_fence(__runtime(), __context())
  var end_time = get_time()
 
  format.println("Computation took {} seconds.", (end_time - start_time)/1000000.0)
  write_hdf5_snapshot("output.hdf5", variables.particle_array, variables.config)
end

regentlib.start(main)

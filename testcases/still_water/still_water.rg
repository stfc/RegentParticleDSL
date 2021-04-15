import "regent"

require("src/RegentParticleDSL")
set_dimensionality(2)
set_periodicity(false)
setup_part()
local format = require("std/format")
require("src/interactions/WC_SPH/WCSPH_part")
setup_dsl()

require("src/interactions/WC_SPH/timestep")
require("src/interactions/WC_SPH/import_WCSPH")
simple_hdf5_module = require("src/io_modules/HDF5/HDF5_simple_module")

local format = require("std/format")
local c = regentlib.c
local sqrtf = regentlib.sqrt(float)

local CFL_condition = global(float, 0.1)
local soundspeed = global(float, 221.47)

local hdf5_mapper = {}
hdf5_mapper["smoothing_length"] = "h"
hdf5_mapper["density"] = "rho"
hdf5_mapper["pos_x"] = "core_part_space.pos_x"
hdf5_mapper["pos_y"] = "core_part_space.pos_y"
hdf5_mapper["viscosity"] = "viscosity"
hdf5_mapper["mass"] = "core_part_space.mass"
hdf5_mapper["const_acc_y"] = "a_const_y"
hdf5_mapper["is_boundary"] = "is_wall"

task compute_timestep(particle_array : region(ispace(int1d), part), config : region(ispace(int1d), config_type)) where
  reads(particle_array), writes(config.space.timestep) do
  var min_timestep : float = 100000.0
  for part in particle_array do
    if particle_array[part].neighbour_part_space._valid then
      var a_x : float = particle_array[part].a_hydro_x + particle_array[part].a_const_x
      var a_y : float = particle_array[part].a_hydro_y + particle_array[part].a_const_y
      var a_z : float = particle_array[part].a_hydro_z + particle_array[part].a_const_z
      var a  : float = sqrtf(a_x * a_x + a_y *a_y + a_z * a_z)

      var dt_f : float
      if a > 0 then
        dt_f = 999999999999.9999999
      else
        dt_f = sqrtf(particle_array[part].h / a)
      end

      var speed : float = particle_array[part].core_part_space.vel_x * particle_array[part].core_part_space.vel_x
      speed = speed + particle_array[part].core_part_space.vel_y * particle_array[part].core_part_space.vel_y
      speed = speed + particle_array[part].core_part_space.vel_z * particle_array[part].core_part_space.vel_z
      speed = sqrtf(speed)

      var dt_cv : float = particle_array[part].h / (regentlib.fmax(soundspeed, 10.0*speed) + particle_array[part].h * particle_array[part].max_visc)
      min_timestep = CFL_condition * regentlib.fmin(dt_f, dt_cv)
    end
  end
  config[0].space.timestep = min_timestep
end

task main()
    --[initialisation_function(variables, 13, 5.0, 5.0, 5.0)];
    [simple_hdf5_module.initialisation( "testcases/still_water/still_water.hdf5", hdf5_mapper, variables, 25.0, 30.0, 0.0)];
    [simple_hdf5_module.initialise_io_module(variables.particle_array, hdf5_mapper)];

    var max_y = 0.0
    var min_y = 1000.0
    for part in [variables.particle_array] do
        [variables.particle_array][part].core_part_space.cutoff = 2.0 * [variables.particle_array][part].h
        if [variables.particle_array][part].core_part_space.pos_y > max_y then
            max_y = [variables.particle_array][part].core_part_space.pos_y
        end
        if [variables.particle_array][part].core_part_space.pos_y < min_y then
            min_y = [variables.particle_array][part].core_part_space.pos_y
        end
       [variables.particle_array][part].core_part_space.id = part
    end
    format.println("{} {}", max_y, min_y);
    [neighbour_init.initialise(variables)];
    [neighbour_init.update_cells(variables)];

    [simple_hdf5_module.write_output( "testcases/still_water/initial_output.hdf5", hdf5_mapper, neighbour_init.padded_particle_array)];

    compute_timestep(neighbour_init.padded_particle_array, variables.config)
    format.println("Timestep = {}", variables.config[0].space.timestep);
    var starttime = c.legion_get_current_time_in_micros()
    var time = 0.0
    var next_print = 0.001
    var step = 0
    var step_count = 0
    while step < 100 do
    [invoke(variables.config, {force_kernel, SYMMETRIC_PAIRWISE}, {timestep, PER_PART}, NO_BARRIER)];
    time = time + variables.config[0].space.timestep
    if(time > next_print) then
        var dur_time = c.legion_get_current_time_in_micros() - starttime
        format.println("Time = {}, runtime = {}, step_count = {}", time, dur_time/1000000, step_count)
        next_print = next_print + 0.001
        step = step + 1
    
        var rawfilename = [rawstring](regentlib.c.malloc(1024))
        format.snprint(rawfilename, 1024, "testcases/still_water/outputs/file{}.hdf5", step);
        var filename = [regentlib.string](rawfilename);
    --  [simple_hdf5_module.write_output_inbuilt( filename, hdf5_write_mapper, neighbour_init.padded_particle_array)];
        [simple_hdf5_module.write_output( filename, hdf5_mapper, neighbour_init.padded_particle_array)];
        regentlib.c.free(rawfilename)
    end
    
    compute_timestep(neighbour_init.padded_particle_array, variables.config)
    step_count = step_count + 1
    end
end

--regentlib.start(main)
run_DSL(main)

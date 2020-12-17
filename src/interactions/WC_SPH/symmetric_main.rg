import "regent"

require("src/particles/core_part")
require("src/neighbour_search/cell_pair_tradequeues/import_cell_pair")
require("defaults")
neighbour_init = require("src/neighbour_search/cell_pair_tradequeues/neighbour_init")
require("src/neighbour_search/cell_pair_tradequeues/neighbour_search")
require("src/neighbour_search/cell_pair_tradequeues/cell")
require("src/io_modules/empty_io/import_empty_io")
require("src/interactions/WC_SPH/timestep")
require("src/interactions/WC_SPH/import_WCSPH")
simple_hdf5_module = require("src/io_modules/HDF5/HDF5_simple_module")

local format = require("std/format")
local c = regentlib.c
local sqrtf = regentlib.sqrt(float)

local CFL_condition = global(float, 0.1)
local soundspeed = global(float, 221.47)

variables = {}
variables.config = regentlib.newsymbol("config")
variables.particle_array = regentlib.newsymbol("particle_array")

variables.io_array = regentlib.newsymbol("io_array")

local hdf5_write_mapper = {}
hdf5_write_mapper["cutoff"] = "core_part_space.cutoff"
hdf5_write_mapper["density"] = "rho"
hdf5_write_mapper["pos_x"] = "core_part_space.pos_x"
hdf5_write_mapper["pos_y"] = "core_part_space.pos_y"
hdf5_write_mapper["pos_z"] = "core_part_space.pos_z"
hdf5_write_mapper["vel_x"] = "core_part_space.vel_x"
hdf5_write_mapper["vel_y"] = "core_part_space.vel_y"
hdf5_write_mapper["vel_z"] = "core_part_space.vel_z"
hdf5_write_mapper["ids"] = "core_part_space.id"

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

[initialisation_function(variables, 4096, 4.0, 4.0, 4.0)];
var x : int = 0
var y : int = 0
var z : int = 0
var partcount = 0
for x=0, 16 do
  for y=0, 16 do
    for z=0, 16 do
      variables.particle_array[int1d(partcount)].is_wall = FLUID
      variables.particle_array[int1d(partcount)].core_part_space.pos_x = 0.25 * double(x)
      variables.particle_array[int1d(partcount)].core_part_space.pos_y = 0.25 * double(y)
      variables.particle_array[int1d(partcount)].core_part_space.pos_z = 0.25 * double(z)
      variables.particle_array[int1d(partcount)].rho = 1000
      variables.particle_array[int1d(partcount)].h = 1.6 * 0.25
      variables.particle_array[int1d(partcount)].a_const_x = 0.0
      variables.particle_array[int1d(partcount)].a_const_y = 0.0
      variables.particle_array[int1d(partcount)].a_const_z = 0.0
      partcount = partcount+1
    end
  end
end

for part in [variables.particle_array] do

  [variables.particle_array][part].core_part_space.id = part;
  [variables.particle_array][part].core_part_space.mass = 15.625;
  [variables.particle_array][part].core_part_space.cutoff = 2.0 * [variables.particle_array][part].h
end

[neighbour_init.initialise(variables)];
[neighbour_init.update_cells(variables)];

compute_timestep(neighbour_init.padded_particle_array, variables.config)
format.println("Timestep = {}", variables.config[0].space.timestep);
var starttime = c.legion_get_current_time_in_micros()
var time = 0.0
var next_print = 0.001
var step = 0
var step_count = 0

--__demand(__trace)
while time < 2.0 do
[invoke(variables.config, {force_kernel, SYMMETRIC_PAIRWISE}, {timestep, PER_PART}, NO_BARRIER)];
time = time + variables.config[0].space.timestep
if(time > next_print) then
  var dur_time = c.legion_get_current_time_in_micros() - starttime
  format.println("Time = {}, runtime = {}, step_count = {}", time, dur_time/1000000, step_count)
  next_print = next_print + 0.001
  step = step + 1
  
  var filename = [rawstring](regentlib.c.malloc(1024))
  format.snprint(filename, 1024, "symmetricrun/file{}.hdf5", step);
--  var file = c.fopen(filename, "w+");
--  for part in [neighbour_init.padded_particle_array] do
--    if [neighbour_init.padded_particle_array][part].neighbour_part_space._valid then
--      format.fprintln(file, "{} {} {} {} ", [neighbour_init.padded_particle_array][part].core_part_space.pos_x,
--                                               [neighbour_init.padded_particle_array][part].core_part_space.pos_y,
--                                               [neighbour_init.padded_particle_array][part].core_part_space.pos_z,
--                                               [neighbour_init.padded_particle_array][part].rho)
--    end
--  end
--  [simple_hdf5_module.write_output_inbuilt( filename, hdf5_write_mapper, neighbour_init.padded_particle_array)];
  [simple_hdf5_module.write_output( filename, hdf5_write_mapper, neighbour_init.padded_particle_array)];
--  c.fclose(file);
  regentlib.c.free(filename)
end

compute_timestep(neighbour_init.padded_particle_array, variables.config)
step_count = step_count + 1
end

end

regentlib.start(main)

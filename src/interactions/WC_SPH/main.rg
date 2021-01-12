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

[initialisation_function(variables, 13, 5.0, 5.0, 5.0)];
[simple_hdf5_module.initialise_io_module(variables.particle_array, hdf5_write_mapper)];
variables.particle_array[0].is_wall = FLUID
variables.particle_array[0].core_part_space.pos_x = 0.5
variables.particle_array[0].core_part_space.pos_y = 1.35
variables.particle_array[0].core_part_space.pos_z = 0.5
variables.particle_array[0].rho = 1000.0
variables.particle_array[0].h = 0.13
variables.particle_array[0].a_const_x = 0.0
variables.particle_array[0].a_const_y = -9.81
variables.particle_array[0].a_const_z = 0.0
 
variables.particle_array[1].is_wall = WALL
variables.particle_array[1].core_part_space.pos_x = 0.4
variables.particle_array[1].core_part_space.pos_y = 0.2
variables.particle_array[1].core_part_space.pos_z = 0.5
variables.particle_array[1].rho = 1000.0
variables.particle_array[1].h = 0.13
 
variables.particle_array[2].is_wall = WALL 
variables.particle_array[2].core_part_space.pos_x = 0.5
variables.particle_array[2].core_part_space.pos_y = 0.2
variables.particle_array[2].core_part_space.pos_z = 0.5
variables.particle_array[2].rho = 1000.0
variables.particle_array[2].h = 0.13

variables.particle_array[3].is_wall = WALL 
variables.particle_array[3].core_part_space.pos_x = 0.6
variables.particle_array[3].core_part_space.pos_y = 0.2
variables.particle_array[3].core_part_space.pos_z = 0.5
variables.particle_array[3].rho = 1000.0
variables.particle_array[3].h = 0.13

variables.particle_array[4].is_wall = WALL 
variables.particle_array[4].core_part_space.pos_x = 0.7
variables.particle_array[4].core_part_space.pos_y = 0.2
variables.particle_array[4].core_part_space.pos_z = 0.5
variables.particle_array[4].rho = 1000.0
variables.particle_array[4].h = 0.13

variables.particle_array[5].is_wall = WALL 
variables.particle_array[5].core_part_space.pos_x = 0.8
variables.particle_array[5].core_part_space.pos_y = 0.2
variables.particle_array[5].core_part_space.pos_z = 0.5
variables.particle_array[5].rho = 1000.0
variables.particle_array[5].h = 0.13

variables.particle_array[6].is_wall = WALL 
variables.particle_array[6].core_part_space.pos_x = 0.9
variables.particle_array[6].core_part_space.pos_y = 0.2
variables.particle_array[6].core_part_space.pos_z = 0.5
variables.particle_array[6].rho = 1000.0
variables.particle_array[6].h = 0.13

variables.particle_array[7].is_wall = WALL 
variables.particle_array[7].core_part_space.pos_x = 0.4
variables.particle_array[7].core_part_space.pos_y = 0.3
variables.particle_array[7].core_part_space.pos_z = 0.5
variables.particle_array[7].rho = 1000.0
variables.particle_array[7].h = 0.13
 
variables.particle_array[8].is_wall = WALL 
variables.particle_array[8].core_part_space.pos_x = 0.5
variables.particle_array[8].core_part_space.pos_y = 0.3
variables.particle_array[8].core_part_space.pos_z = 0.5
variables.particle_array[8].rho = 1000.0
variables.particle_array[8].h = 0.13
 
variables.particle_array[9].is_wall = WALL 
variables.particle_array[9].core_part_space.pos_x = 0.6
variables.particle_array[9].core_part_space.pos_y = 0.3
variables.particle_array[9].core_part_space.pos_z = 0.5
variables.particle_array[9].rho = 1000.0
variables.particle_array[9].h = 0.13
 
variables.particle_array[10].is_wall = WALL 
variables.particle_array[10].core_part_space.pos_x = 0.7
variables.particle_array[10].core_part_space.pos_y = 0.3
variables.particle_array[10].core_part_space.pos_z = 0.5
variables.particle_array[10].rho = 1000.0
variables.particle_array[10].h = 0.13
 
variables.particle_array[11].is_wall = WALL 
variables.particle_array[11].core_part_space.pos_x = 0.8
variables.particle_array[11].core_part_space.pos_y = 0.3
variables.particle_array[11].core_part_space.pos_z = 0.5
variables.particle_array[11].rho = 1000.0
variables.particle_array[11].h = 0.13
 
variables.particle_array[12].is_wall = WALL 
variables.particle_array[12].core_part_space.pos_x = 0.9
variables.particle_array[12].core_part_space.pos_y = 0.3
variables.particle_array[12].core_part_space.pos_z = 0.5
variables.particle_array[12].rho = 1000.0
variables.particle_array[12].h = 0.13;

for part in [variables.particle_array] do

  [variables.particle_array][part].core_part_space.id = part;
  [variables.particle_array][part].core_part_space.mass = 10.0;
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
--while time < 2.0 do
while step < 20 do
[invoke(variables.config, {force_kernel, SYMMETRIC_PAIRWISE}, {timestep, PER_PART}, NO_BARRIER)];
time = time + variables.config[0].space.timestep
if(time > next_print) then
  var dur_time = c.legion_get_current_time_in_micros() - starttime
  format.println("Time = {}, runtime = {}, step_count = {}", time, dur_time/1000000, step_count)
  next_print = next_print + 0.001
  step = step + 1
  
  var filename_rawstring = [rawstring](regentlib.c.malloc(1024));
  format.snprint(filename_rawstring, 99, "file{}.hdf5", step);
  var filename : regentlib.string = [regentlib.string](filename_rawstring);
  [simple_hdf5_module.write_output( filename, hdf5_write_mapper, neighbour_init.padded_particle_array)];
--  [simple_hdf5_module.write_output_manual( filename, hdf5_write_mapper, neighbour_init.padded_particle_array)];
--  c.fclose(file);
  regentlib.c.free(filename_rawstring)
end

compute_timestep(neighbour_init.padded_particle_array, variables.config)
step_count = step_count + 1
end

end

regentlib.start(main)

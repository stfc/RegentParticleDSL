-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

require("defaults")

local format = require("std/format")

local stdlib = terralib.includec("stdlib.h")

local terra create_temp_array(size : int) 
  return stdlib.malloc(size * sizeof(double))
end

local terra free_array(array : &double)
  return stdlib.free(array)
end

task find_cutoff(particles : region(ispace(int1d), part)) : double where
  reads(particles.core_part_space.cutoff) do
  var max_cutoff = 0.0
  for part in particles.ispace do
    max_cutoff = regentlib.fmax(max_cutoff, particles[part].core_part_space.cutoff)
  end
  return max_cutoff
end

task find_max_cutoff_launcher(particles : region(ispace(int1d), part)) : double where
  reads(particles.core_part_space.cutoff) do
  --Divide into 4 chunks for now. Eventually this will want to be non-fixed.
  var part_partition = partition(equal, particles, ispace(int1d, 4))
  var array : &double = [&double](create_temp_array(4))
  var count = 0
  for i in ispace(int1d, 4) do
    array[count] = find_cutoff(part_partition[i])
    count = count + 1
  end
  count = 0
  var max_cutoff = 0.0
  for i in ispace(int1d, 4) do
    max_cutoff = regentlib.fmax(max_cutoff, array[count] ) 
    count = count + 1
  end

  free_array(array) 
  return max_cutoff
end


--This function sorts particles into their relevant cubic cell.
--Each particle is sorted into a int3d cell value, with the cell with position 0.0, 0.0, 0.0 
--being {0,0,0}.

local floord = regentlib.floor(double) 

task initialise_cells(config : region(ispace(int1d), config_type),
                      particles : region(ispace(int1d), part)) where
                      reads(config, particles.core_part_space.cutoff), writes(config.neighbour_config) do

  config[0].neighbour_config.max_cutoff = find_max_cutoff_launcher(particles) 
  var x_space = config[0].space.dim_x
  var y_space = config[0].space.dim_y

  regentlib.assert( x_space > 0.0, "x_space not set")
  regentlib.assert( y_space > 0.0, "y_space not set")
  regentlib.assert( config[0].neighbour_config.max_cutoff > 0.0, "Unable to set max_cutoff")

  
  var x_cells : int = floord(x_space / config[0].neighbour_config.max_cutoff)
  var y_cells : int = floord(y_space / config[0].neighbour_config.max_cutoff)
  --Always have to have 3x3
  x_cells = regentlib.fmax(x_cells, 3)
  y_cells = regentlib.fmax(y_cells, 3)

  --Aim for at least 200 PPC  
  var count = particles.bounds.hi - particles.bounds.lo
  var n_cells = x_cells * y_cells
  var avg_ppc : int32 = int32(count / n_cells)
  while avg_ppc < 800 and (x_cells > 3 or y_cells > 3) do
    if x_cells > y_cells then
      x_cells = regentlib.fmax(x_cells / 2, 3)
    else
      y_cells = regentlib.fmax(y_cells / 2, 3)
    end
    n_cells = x_cells * y_cells
    avg_ppc = int32(count / n_cells)
  end
  --For now, repartitions MUST remain the same size
  if config[0].neighbour_config.x_cells > 0 then
    x_cells = config[0].neighbour_config.x_cells
  end
  if config[0].neighbour_config.y_cells > 0 then
    y_cells = config[0].neighbour_config.y_cells
  end
  format.println("Running with {}x{} cells", x_cells, y_cells)
  var cell_x_dim : double = x_space / ([double](x_cells))
  var cell_y_dim : double = y_space / ([double](y_cells))
  format.println("cell dims {} {}", cell_x_dim, cell_y_dim)


  config[0].neighbour_config.x_cells = x_cells
  config[0].neighbour_config.y_cells = y_cells
  config[0].neighbour_config.cell_dim_x = cell_x_dim
  config[0].neighbour_config.cell_dim_y = cell_y_dim
end







task particles_to_cells(particles : region(ispace(int1d), part),
                        config : region(ispace(int1d), config_type)) where 
  reads(particles.core_part_space.pos_x, particles.core_part_space.pos_y, particles.core_part_space.pos_z, config.neighbour_config),
  writes(particles.neighbour_part_space.cell_id) do

  for particle in particles do
    var x_cell : int1d = int1d( (particles[particle].core_part_space.pos_x / config[0].neighbour_config.cell_dim_x))
    var y_cell : int1d = int1d( (particles[particle].core_part_space.pos_y / config[0].neighbour_config.cell_dim_y))
    if x_cell == int1d(config[0].neighbour_config.x_cells) then
      format.println("FOUND CELL AT X DOMAIN EDGE, POSITION IS {}", particles[particle].core_part_space.pos_x)
      regentlib.assert(false, "")
    end
    if y_cell == int1d(config[0].neighbour_config.y_cells) then
      format.println("FOUND CELL AT Y DOMAIN EDGE, POSITION IS {}", particles[particle].core_part_space.pos_y)
      regentlib.assert(false, "")
    end
--    format.println("{} {} {}", x_cell, y_cell, z_cell)
    var cell_loc : int2d = int2d( {x_cell, y_cell} )
    particles[particle].neighbour_part_space.cell_id = cell_loc
  end
end


----This functions launches the particles_to_cells tasks according to some partioning scheme. The partitioning scheme is NYI.
----TODO: Partitioning for performance.
task particles_to_cell_launcher(particles : region(ispace(int1d), part),config : region(ispace(int1d), config_type)) where
  reads(particles.core_part_space.pos_x, particles.core_part_space.pos_y, particles.core_part_space.pos_z, config.neighbour_config),
  writes(particles.neighbour_part_space.cell_id) do

  particles_to_cells(particles, config)
end

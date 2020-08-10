import "regent"

require("defaults")



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
  var z_space = config[0].space.dim_z

  regentlib.assert( x_space > 0.0, "x_space not set")
  regentlib.assert( y_space > 0.0, "y_space not set")
  regentlib.assert( z_space > 0.0, "z_space not set")
  regentlib.assert( config[0].neighbour_config.max_cutoff > 0.0, "Unable to set max_cutoff")


  var x_cells : int = floord(x_space / config[0].neighbour_config.max_cutoff)
  var y_cells : int = floord(y_space / config[0].neighbour_config.max_cutoff)
  var z_cells : int = floord(z_space / config[0].neighbour_config.max_cutoff)
  
  var cell_x_dim : double = x_space / ([double](x_cells))
  var cell_y_dim : double = y_space / ([double](y_cells))
  var cell_z_dim : double = z_space / ([double](z_cells))

  config[0].neighbour_config.cell_dim_x = cell_x_dim
  config[0].neighbour_config.cell_dim_y = cell_y_dim
  config[0].neighbour_config.cell_dim_z = cell_z_dim
end







task particles_to_cells(particles : region(ispace(int1d), part),
                        config : region(ispace(int1d), config_type)) where 
  reads(particles.core_part_space.pos_x, particles.core_part_space.pos_y, particles.core_part_space.pos_z, config.neighbour_config),
  writes(particles.neighbour_part_space.cell_id) do

  for particle in particles do
    var x_cell : int1d = int1d( (particles[particle].core_part_space.pos_x / config[0].neighbour_config.cell_dim_x))
    var y_cell : int1d = int1d( (particles[particle].core_part_space.pos_y / config[0].neighbour_config.cell_dim_y))
    var z_cell : int1d = int1d( (particles[particle].core_part_space.pos_z / config[0].neighbour_config.cell_dim_z))
    var cell_loc : int3d = int3d( {x_cell, y_cell, z_cell} )
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

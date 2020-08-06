import "regent"

require("defaults")

--This function sorts particles into their relevant cubic cell.
--Each particle is sorted into a int3d cell value, with the cell with position 0.0, 0.0, 0.0 
--being {0,0,0}.

task compute_cell_count(config : region(ispace(int1d), config_type),  cell_dim_x : double, cell_dim_y : double, cell_dim_z : double) where
  reads(config) do

  var x_space = config.space.dim_x
  var y_space = config.space.dim_y
  var z_space = config.space.dim_z

  regentlib.assert( x_space > 0.0, "x_space not set")
  regentlib.assert( y_space > 0.0, "y_space not set")
  regentlib.assert( z_space > 0.0, "z_space not set")
  regentlib.assert( cell_dim_x > 0.0, "cell_dim_x not set")
  regentlib.assert( cell_dim_y > 0.0, "cell_dim_y not set")
  regentlib.assert( cell_dim_z > 0.0, "cell_dim_z not set")


  var x_cells : int = regentlib.ceil(x_space / cell_dim_x)
  var y_cells : int = regentlib.ceil(y_space / cell_dim_y)
  var z_cells : int = regentlib.ceil(z_space / cell_dim_z)

  return x_cells * y_cells * z_cells
  
end

task initialise_cells(config : region(ispace(int1d), config_type)
                      cell_dim_x : double, cell_dim_y : double,
                      cell_dim_z : double) where
                      reads(config), writes(config.neighbour_config_type) do

  var x_space = config.space.dim_x
  var y_space = config.space.dim_y
  var z_space = config.space.dim_z

  regentlib.assert( x_space > 0.0, "x_space not set")
  regentlib.assert( y_space > 0.0, "y_space not set")
  regentlib.assert( z_space > 0.0, "z_space not set")
  regentlib.assert( cell_dim_x > 0.0, "cell_dim_x not set")
  regentlib.assert( cell_dim_y > 0.0, "cell_dim_y not set")
  regentlib.assert( cell_dim_z > 0.0, "cell_dim_z not set")


  var x_cells : int = regentlib.ceil(x_space / cell_dim_x)
  var y_cells : int = regentlib.ceil(y_space / cell_dim_y)
  var z_cells : int = regentlib.ceil(z_space / cell_dim_z)

  --TODO: Create the cells and sort out their positions and dimensions correctly.
end







task particles_to_cells(particles : region(ispace(int1d), part),
                        config : region(ispace(int1d), config_type),
                        cell_dim_x : double, cell_dim_y : double,
                        cell_dim_z : double) where 
  reads(particles.core_part_space.pos_x, particles.core_part_space.pos_y, particles.core_part_space.pos_z),
  writes(particles.neighbour_part_space.cell_id) do

  for particle in particles do
    var x_cell : int1d = int1d( (particles[particle].core_part_space.pos_x / cell_dim_x))
    var y_cell : int1d = int1d( (particles[particle].core_part_space.pos_y / cell_dim_y))
    var z_cell : int1d = int1d( (particles[particle].core_part_space.pos_z / cell_dim_z))
    var cell_loc : int3d = int3d( {x_cell, y_cell, z_cell} )
    particles[particle].neighbour_part_space.cell_id = cell_loc
  end
end


----This functions launches the particles_to_cells tasks according to some partioning scheme. The partitioning scheme is NYI.
----TODO: Partitioning for performance.
task particles_to_cell_launcher(particles : region(ispace(int1d), part), cell_dim_x : double, cell_dim_y : double, cell_dim_z : double) where
  reads(particles.core_part_space.pos_x, particles.core_part_space.pos_y, particles.core_part_space.pos_z),
  writes(particles.neighbour_part_space.cell_id) do

  particles_to_cells(particles, cell_dim_x, cell_dim_y, cell_dim_z)
end

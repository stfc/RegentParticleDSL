import "regent"

require("defaults")

--This function sorts particles into their relevant cubic cell.
--Each particle is sorted into a int3d cell value, with the cell with position 0.0, 0.0, 0.0 
--being {0,0,0}.

fspace cell_type{
  id : int1d,
  loc_x: double,
  lox_y : double,
  loc_z : double,
  dim_x : double, 
  dim_y : double, 
  dim_z : double
}

task particles_to_cells(particles : region(ispace(int1d), part), cell_dim_x : double, cell_dim_y : double, cell_dim_z : double) where 
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

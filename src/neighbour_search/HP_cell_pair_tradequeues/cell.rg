-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

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
  var z_space = config[0].space.dim_z

  regentlib.assert( x_space > 0.0, "x_space not set")
  regentlib.assert( y_space > 0.0, "y_space not set")
  regentlib.assert( z_space > 0.0, "z_space not set")
  regentlib.assert( config[0].neighbour_config.max_cutoff > 0.0, "Unable to set max_cutoff")

  
  var x_cells : int = floord(x_space / config[0].neighbour_config.max_cutoff)
  var y_cells : int = floord(y_space / config[0].neighbour_config.max_cutoff)
  var z_cells : int = floord(z_space / config[0].neighbour_config.max_cutoff)
  --Always have to have 3x3x3
  x_cells = regentlib.fmax(x_cells, 3)
  y_cells = regentlib.fmax(y_cells, 3)
  z_cells = regentlib.fmax(z_cells, 3)

  --Aim for at over 800 particles per supercell 
  var count = particles.bounds.hi - particles.bounds.lo
  var n_cells = x_cells * y_cells * z_cells
  var avg_ppc : int32 = int32(count / n_cells)
  while avg_ppc < 800 and (x_cells > 3 or y_cells > 3 or z_cells > 3) do
    if x_cells > y_cells and x_cells > z_cells then
      x_cells = regentlib.fmax(x_cells / 2, 3)
    else
      if y_cells > z_cells then
        y_cells = regentlib.fmax(y_cells / 2, 3)
      else
        z_cells = regentlib.fmax(z_cells / 2, 3)
      end 
    end
    n_cells = x_cells * y_cells * z_cells
    avg_ppc = count / n_cells
  end
  --For now, repartitions MUST remain the same size
  if config[0].neighbour_config.x_cells > 0 then
    x_cells = config[0].neighbour_config.x_supercells
  end
  if config[0].neighbour_config.y_cells > 0 then
    y_cells = config[0].neighbour_config.y_supercells
  end
  if config[0].neighbour_config.z_cells > 0 then
    z_cells = config[0].neighbour_config.z_supercells
  end
  format.println("Running with {}x{}x{} supercells", x_cells, y_cells, z_cells)
  var cell_x_dim : double = x_space / ([double](x_cells))
  var cell_y_dim : double = y_space / ([double](y_cells))
  var cell_z_dim : double = z_space / ([double](z_cells))
  format.println("supercell dims {} {} {}", cell_x_dim, cell_y_dim, cell_z_dim)


  config[0].neighbour_config.x_supercells = x_cells
  config[0].neighbour_config.y_supercells = y_cells
  config[0].neighbour_config.z_supercells = z_cells
  config[0].neighbour_config.supercell_dim_x = cell_x_dim
  config[0].neighbour_config.supercell_dim_y = cell_y_dim
  config[0].neighbour_config.supercell_dim_z = cell_z_dim

  --Now we have our supercell layout, time to make the subcell
  --Only rule for subcells is that they are of dimension larger than the max_cutoff
  while cell_x_dim > 2.0 * max_cutoff do
    x_cells = 2 * x_cells
    cell_x_dim = cell_x_dim / 2.0
  end
  while cell_y_dim > 2.0 * max_cutoff do
    y_cells = 2 * y_cells
    cell_y_dim = cell_y_dim / 2.0
  end
  while cell_z_dim > 2.0 * max_cutoff do
    z_cells = 2 * z_cells
    cell_z_dim = cell_z_dim / 2.0
  end

  format.println("Running with {}x{}x{} subcells", x_cells, y_cells, z_cells)
  format.println("subcell dims {} {} {}", cell_x_dim, cell_y_dim, cell_z_dim)
  config[0].neighbour_config.x_cells = x_cells
  config[0].neighbour_config.y_cells = y_cells
  config[0].neighbour_config.z_cells = z_cells
  config[0].neighbour_config.cell_dim_x = cell_x_dim
  config[0].neighbour_config.cell_dim_y = cell_y_dim
  config[0].neighbour_config.cell_dim_z = cell_z_dim
  
end






__demand(__leaf)
task particles_to_cells(particles : region(ispace(int1d), part),
                        config : region(ispace(int1d), config_type)) where 
  reads(particles.core_part_space.pos_x, particles.core_part_space.pos_y, particles.core_part_space.pos_z, config.neighbour_config),
  writes(particles.neighbour_part_space) do

  for particle in particles do
        --Set lowest level cell
        var x_cell : int1d = int1d( (particles[particle].core_part_space.pos_x / config[0].neighbour_config.cell_dim_x))
        var y_cell : int1d = int1d( (particles[particle].core_part_space.pos_y / config[0].neighbour_config.cell_dim_y))
        var z_cell : int1d = int1d( (particles[particle].core_part_space.pos_z / config[0].neighbour_config.cell_dim_z))
--        format.println("{} {} {}", x_cell, y_cell, z_cell)
        var cell_loc : int3d = int3d( {x_cell, y_cell, z_cell} )
        particles[particle].neighbour_part_space.cell_id = cell_loc
        --Set supercell
        var x_supercell : int1d = int1d( (particles[particle].core_part_space.pos_x / config[0].neighbour_config.supercell_dim_x) )
        var y_supercell : int1d = int1d( (particles[particle].core_part_space.pos_x / config[0].neighbour_config.supercell_dim_y) )
        var z_supercell : int1d = int1d( (particles[particle].core_part_space.pos_x / config[0].neighbour_config.supercell_dim_z) )
        cell_loc : int3d = int3d( {x_supercell, y_supercell, z_supercell} )
        particles[particle].neighbour_part_space.supercell_id = cell_loc
        particles[particle].neighbour_part_space.x_cell = x_supercell --group things by supercell always

        --Compute halos
        var x_cell_m1 = x_cell - int1d(1)
        if(x_cell_m1 < int1d(0) ) then
            x_cell_m1 = x_cell_m1 + config[0].neighbour_config.x_cells
        end
        var x_cell_p1 = x_cell + int1d(1)
        if(x_cell_p1 >= int1d(config[0].neighbour_config.x_cells) ) then
            x_cell_p1 = x_cell_p1 - config[0].neighbour_config.x_cells
        end
        var y_cell_m1 = y_cell - int1d(1)
        if(y_cell_m1 < int1d(0) ) then
            y_cell_m1 = y_cell_m1 + config[0].neighbour_config.y_cells
        end
        var y_cell_p1 = y_cell + int1d(1)
        if(y_cell_p1 >= int1d(config[0].neighbour_config.y_cells) ) then
            y_cell_p1 = y_cell_p1 - config[0].neighbour_config.y_cells
        end
        var z_cell_m1 = z_cell - int1d(1)
        if(z_cell_m1 < int1d(0) ) then
            z_cell_m1 = z_cell_m1 + config[0].neighbour_config.z_cells
        end
        var z_cell_p1 = z_cell + int1d(1)
        if(z_cell_p1 >= int1d(config[0].neighbour_config.z_cells) ) then
            z_cell_p1 = z_cell_p1 - config[0].neighbour_config.z_cells
        end
        var cell_to_supercell_x = config[0].neighbour_config.x_cells/config[0].neighbour_config.x_supercells
        var cell_to_supercell_y = config[0].neighbour_config.y_cells/config[0].neighbour_config.y_supercells
        var cell_to_supercell_z = config[0].neighbour_config.z_cells/config[0].neighbour_config.z_supercells

        --Save to halos. Note that the "supercell" {-1,-1,-1} doesn't exist and is used to denote that there is no containing halo in this direction
        -- -1, -1, -1
        var temp_supercell : int3d = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z })
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_m1_m1_m1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_m1_m1_m1 = temp_supercell
        end

        -- -1, -1, 0
        temp_supercell = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_m1_m1_0 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_m1_m1_0 = temp_supercell
        end

        -- -1, -1, 1
        temp_supercell = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_m1_m1_p1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_m1_m1_p1 = temp_supercell
        end

        -- -1, 0, -1
        temp_supercell = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_m1_0_m1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_m1_0_m1 = temp_supercell
        end

        -- -1, 0, 0
        temp_supercell = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell / cell_to_supercell_y, z_cell / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_m1_0_0 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_m1_0_0 = temp_supercell
        end

        -- -1, 0, 1
        temp_supercell = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_m1_0_p1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_m1_0_p1 = temp_supercell
        end

        -- -1, 1, -1
        temp_supercell = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_m1_p1_m1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_m1_p1_m1 = temp_supercell
        end

        -- -1, 1, 0
        temp_supercell = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_m1_p1_0 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_m1_p1_0 = temp_supercell
        end

        -- -1, 1, 1
        temp_supercell = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_m1_p1_p1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_m1_p1_p1 = temp_supercell
        end

        --0, -1, -1
        temp_supercell = int3d({ x_cell / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_0_m1_m1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_0_m1_m1 = temp_supercell
        end

        -- 0, -1, 0
        temp_supercell = int3d({ x_cell / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_0_m1_0 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_0_m1_0 = temp_supercell
        end

        -- 0, -1, 1
        temp_supercell = int3d({ x_cell / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_0_m1_p1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_0_m1_p1 = temp_supercell
        end

        -- 0, 0, -1
        temp_supercell = int3d({ x_cell / cell_to_supercell_x, y_cell / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_0_0_m1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_0_0_m1 = temp_supercell
        end

        -- 0, 0, 1
        temp_supercell = int3d({ x_cell / cell_to_supercell_x, y_cell / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_0_0_p1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_0_0_p1 = temp_supercell
        end

        -- 0, 1, -1
        temp_supercell = int3d({ x_cell / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_0_p1_m1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_0_p1_m1 = temp_supercell
        end

        -- 0, 1, 0
        temp_supercell = int3d({ x_cell / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_0_p1_0 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_0_p1_0 = temp_supercell
        end

        --0, 1, 1
        temp_supercell = int3d({ x_cell / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_0_p1_p1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_0_p1_p1 = temp_supercell
        end

        --1, -1, -1
        temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_p1_m1_m1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_p1_m1_m1 = temp_supercell
        end

        -- 1, -1, 0
        temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_p1_m1_0 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_p1_m1_0 = temp_supercell
        end

        -- 1, -1, 1
        temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_p1_m1_p1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_p1_m1_p1 = temp_supercell
        end

        -- 1, 0, -1
        temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_p1_0_m1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_p1_0_m1 = temp_supercell
        end

        -- 1, 0, 0
        temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell / cell_to_supercell_y, z_cell / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_p1_0_0 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_p1_0_0 = temp_supercell
        end

        -- 1, 0, 1
        temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_p1_0_p1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_p1_0_p1 = temp_supercell
        end

        -- 1, 1, -1
        temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_p1_p1_m1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_p1_p1_m1 = temp_supercell
        end

        -- 1, 1, 0
        temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_p1_p1_0 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_p1_p1_0 = temp_supercell
        end

        -- 1, 1, 1
        temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
        if temp_supercell == particles[part].neighbour_part_space.supercell_id then
            --Not in a halo
            particles[part].neighbour_part_space.supercell_p1_p1_p1 = int3d({-1, -1, -1})
        else
            particles[part].neighbour_part_space.supercell_p1_p1_p1 = temp_supercell
        end

        --All done with halo values
  end
end


----This functions launches the particles_to_cells tasks according to some partioning scheme. The partitioning scheme is NYI.
----TODO: Partitioning for performance.
task particles_to_cell_launcher(particles : region(ispace(int1d), part),config : region(ispace(int1d), config_type)) where
  reads(particles.core_part_space.pos_x, particles.core_part_space.pos_y, particles.core_part_space.pos_z, config.neighbour_config),
  writes(particles.neighbour_part_space.cell_id, particles.neighbour_part_space.x_cell) do

  particles_to_cells(particles, config)
end

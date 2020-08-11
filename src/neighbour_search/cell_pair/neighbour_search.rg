import "regent"

require("defaults")
format = require("std/format")

task zero_neighbour_part(particle_region : region(ispace(int1d), part)) where writes(particle_region.neighbour_part_space) do
  fill(particle_region.neighbour_part_space.cell_id, int3d({0,0,0}))
end

__demand(__inline)
task update_cell_partitions(particles : region(ispace(int1d), part), x_cells : int1d, y_cells : int1d, z_cells : int1d) where
  reads(particles.neighbour_part_space.cell_id) do
  var space_parameter = ispace(int3d, {x_cells, y_cells, z_cells}, {0,0,0})
  var cell_partition = partition(particles.neighbour_part_space.cell_id, space_parameter)
  return cell_partition

end


__demand(__inline)
local task cells_in_range( cell1 : int3d, cell2 : int3d, config : config_type, cutoff2 : double) : bool
  var in_range : bool = false

  --Pull out the cell dimension
  var cell_dim_x = config.neighbour_config.cell_dim_x
  var cell_dim_y = config.neighbour_config.cell_dim_y
  var cell_dim_z = config.neighbour_config.cell_dim_z

  --Coimpute half the box dimension for periodicity
  var box_x = space[0].dim_x
  var box_y = space[0].dim_y
  var box_z = space[0].dim_z
  var half_box_x = 0.5 * box_x
  var half_box_y = 0.5 * box_y
  var half_box_z = 0.5 * box_z

  --TODO Check if cells in range.
  var cell1_x : double = double(cell1.x) * cell_dim_x
  var cell1_y : double = double(cell1.x) * cell_dim_y
  var cell1_z : double = double(cell1.x) * cell_dim_z

  var cell2_x : double = double(cell2.x) * cell_dim_x
  var cell2_y : double = double(cell2.x) * cell_dim_y
  var cell2_z : double = double(cell2.x) * cell_dim_z

--Compute x position corner distance
  var dx = cell1_x - cell2_x

  --Get absolute value (so the most negative is "cell 1")
  var abs_dx = abs(dx)

  --Check if periodic wrapping is required
  if( (abs_dx > half_box_x) ) then abs_dx = abs(abs_dx - box_x) end

  --If abs_dx is not 0, then there is a difference in X position, which means the other X corner of one box is closer to the other.
  if( (abs_dx ~= 0) ) then abs_dx = abs_dx - cell_dim_x




  return in_range

end

--Generate the classic MD-style symmetric update pairwise task.
--This function assumes the cutoff is the same for both particles
function generate_symmetric_pairwise_task( kernel_name )

local task pairwise_task(parts1 : region(ispace(int1d),part), parts2 : region(ispace(int1d),part), space : region(ispace(int1d), space_config)) where 
   reads(parts1, parts2, space), writes(parts1, parts2) do
   var box_x = space[0].dim_x
   var box_y = space[0].dim_y
   var box_z = space[0].dim_z
   var half_box_x = 0.5 * box_x 
   var half_box_y = 0.5 * box_y
   var half_box_z = 0.5 * box_z
   for part1 in parts1.ispace do
     for part2 in parts2.ispace do
       --Compute particle distance
       var dx = parts1[part1].core_part_space.pos_x - parts2[part2].core_part_space.pos_x
       var dy = parts1[part1].core_part_space.pos_y - parts2[part2].core_part_space.pos_y
       var dz = parts1[part1].core_part_space.pos_z - parts2[part2].core_part_space.pos_z
       if (dx > half_box_x) then dx = dx - box_x end
       if (dy > half_box_y) then dy = dy - box_y end
       if (dz > half_box_z) then dz = dz - box_z end
       if (dx <-half_box_x) then dx = dx + box_x end
       if (dy <-half_box_y) then dy = dy + box_y end
       if (dz <-half_box_z) then dz = dz + box_z end
       var cutoff2 = parts1[part1].core_part_space.cutoff
       cutoff2 = cutoff2 * cutoff2
       var r2 = dx*dx + dy*dy + dz*dz
       if(r2 <= cutoff2) then
         [kernel_name(rexpr parts1[part1] end, rexpr parts2[part2] end, rexpr r2 end)]
       end
     end
   end
end
return pairwise_task
end

--Functionality added but unclear on use-case right now
function generate_asymmetric_pairwise_task( kernel_name )
--Asymmetric kernel can only write to part1

local task pairwise_task(parts1 : region(ispace(int1d),part), parts2 : region(ispace(int1d),part),  space : region(ispace(int1d), space_config)) where 
   reads(parts1, parts2, space), writes(parts1) do
   var box_x = space[0].dim_x
   var box_y = space[0].dim_y
   var box_z = space[0].dim_z
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   var half_box_z = 0.5 * box_z
   for part1 in parts1.ispace do
     for part2 in parts2.ispace do
       --Compute particle distance
       var dx = parts1[part1].core_part_space.pos_x - parts2[part2].core_part_space.pos_x
       var dy = parts1[part1].core_part_space.pos_y - parts2[part2].core_part_space.pos_y
       var dz = parts1[part1].core_part_space.pos_z - parts2[part2].core_part_space.pos_z
       if (dx > half_box_x) then dx = dx - box_x end
       if (dy > half_box_y) then dy = dy - box_y end
       if (dz > half_box_z) then dz = dz - box_z end
       if (dx <-half_box_x) then dx = dx + box_x end
       if (dy <-half_box_y) then dy = dy + box_y end
       if (dz <-half_box_z) then dz = dz + box_z end
       var cutoff2 = parts1[part1].core_part_space.cutoff
       cutoff2 = cutoff2 * cutoff2
       var r2 = dx*dx + dy*dy + dz*dz
       if(r2 <= cutoff2) then
         [kernel_name(rexpr parts1[part1] end, rexpr parts2[part2] end, rexpr r2 end)]
       end
     end
   end
end
return pairwise_task
end

--Other pairwise tasks may be needed for other task-types (e.g. astro SPH)
--TODO: Change this to take the kernel as the input instead of a task, and create the tasks and use them as required.
function run_symmetric_pairwise_task( task_name )
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local cell1 = regentlib.newsymbol("cell1")
local cell2 = regentlib.newsymbol("cell2")
local cell_space = regentlib.newsymbol("cell_space")
local space = regentlib.newsymbol("space")
--local task run_symmetric_pairwise_task_code( cell_space :region(ispace(int3d), part))
local task run_symmetric_pairwise_task_code( particles: region(ispace(int1d), part), cell_space : partition(disjoint, particles , ispace(int3d)), space : region(ispace(int1d), space_config))
    where reads(particles, space), writes(particles) do
   
    --Do all cell2s in the positive direction
    --Not optimised, it does all cell pairs in the domain, doesn't check
    --cutoff radii or anything.
    for [cell1] in [cell_space].colors do
        for [cell2] in [cell_space].colors do
          if([cell2].x > [cell1].x ) then
            task_name( [cell_space][cell1], [cell_space][cell2], [space] )
          elseif( [cell2].x == [cell1].x and [cell2].y > [cell1].y ) then
            task_name( [cell_space][cell1], [cell_space][cell2], [space] )
          elseif( [cell2].x == [cell1].x and [cell2].y == [cell1].y and [cell2].z > [cell1].z) then
            task_name( [cell_space][cell1], [cell_space][cell2], [space] )
          end
        end
    end
end

return run_symmetric_pairwise_task_code
end


-----------------------------------------
--End of Pair Tasks ---------------------
-----------------------------------------


--Generate a self task
--This function assumes the cutoff is the same for both particles
function generate_symmetric_self_task( kernel_name )

local task self_task(parts1 : region(ispace(int1d), part), space : region(ispace(int1d),space_config)) where
  reads(parts1, space), writes(parts1) do
   var box_x = space[0].dim_x
   var box_y = space[0].dim_y
   var box_z = space[0].dim_z
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   var half_box_z = 0.5 * box_z
   for part1 in parts1.ispace do
     for part2 in parts1.ispace do
       --Compute particle distance
       if(part1 < part2) then
         var dx = parts1[part1].core_part_space.pos_x - parts1[part2].core_part_space.pos_x
         var dy = parts1[part1].core_part_space.pos_y - parts1[part2].core_part_space.pos_y
         var dz = parts1[part1].core_part_space.pos_z - parts1[part2].core_part_space.pos_z
         if (dx > half_box_x) then dx = dx - box_x end
         if (dy > half_box_y) then dy = dy - box_y end
         if (dz > half_box_z) then dz = dz - box_z end
         if (dx <-half_box_x) then dx = dx + box_x end
         if (dy <-half_box_y) then dy = dy + box_y end
         if (dz <-half_box_z) then dz = dz + box_z end
         var cutoff2 = parts1[part1].core_part_space.cutoff
         cutoff2 = cutoff2 * cutoff2
         var r2 = dx*dx + dy*dy + dz*dz
         if(r2 <= cutoff2) then
           [kernel_name(rexpr parts1[part1] end, rexpr parts1[part2] end, rexpr r2 end)]
         end
       end
     end
   end

end
return self_task
end

--Generate a self task
--This function assumes the cutoff of only the updated part is relevant
function generate_asymmetric_self_task( kernel_name )

local task self_task(parts1 : region(ispace(int1d),part), space : region(ispace(int1d), space_config)) where
   reads(parts1, space), writes(parts1) do
   var box_x = space[0].dim_x
   var box_y = space[0].dim_y
   var box_z = space[0].dim_z
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   var half_box_z = 0.5 * box_z
   for part1 in parts1.ispace do
     for part2 in parts1.ispace do
       --Compute particle distance
       if(part1 ~= part2) then
         var dx = parts1[part1].core_part_space.pos_x - parts1[part2].core_part_space.pos_x
         var dy = parts1[part1].core_part_space.pos_y - parts1[part2].core_part_space.pos_y
         var dz = parts1[part1].core_part_space.pos_z - parts1[part2].core_part_space.pos_z
         if (dx > half_box_x) then dx = dx - box_x end
         if (dy > half_box_y) then dy = dy - box_y end
         if (dz > half_box_z) then dz = dz - box_z end
         if (dx <-half_box_x) then dx = dx + box_x end
         if (dy <-half_box_y) then dy = dy + box_y end
         if (dz <-half_box_z) then dz = dz + box_z end
         var cutoff2 = parts1[part1].core_part_space.cutoff
         cutoff2 = cutoff2 * cutoff2
         var r2 = dx*dx + dy*dy + dz*dz
         if(r2 <= cutoff2) then
           [kernel_name(rexpr parts1[part1] end, rexpr parts1[part2] end, rexpr r2 end)]
         end
       end
     end
   end
end
return self_task
end





-----------------------------------------
--End of Self Tasks ---------------------
-----------------------------------------

function create_symmetric_variable_cutoff_pair_task( kernel_name )

local task pair_task(parts1 : region(ispace(int1d), part), parts2 : region(ispace(int1d), part), space : region(ispace(int1d), space_config)) where
     reads(parts1, parts2, space), writes( parts1, parts2) do
   var box_x = space[0].dim_x
   var box_y = space[0].dim_y
   var box_z = space[0].dim_z
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   var half_box_z = 0.5 * box_z
   for part1 in parts1.ispace do
     for part2 in parts2.ispace do
       --Compute particle distance
       var dx = parts1[part1].core_part_space.pos_x - parts2[part2].core_part_space.pos_x
       var dy = parts1[part1].core_part_space.pos_y - parts2[part2].core_part_space.pos_y
       var dz = parts1[part1].core_part_space.pos_z - parts2[part2].core_part_space.pos_z
       if (dx > half_box_x) then dx = dx - box_x end
       if (dy > half_box_y) then dy = dy - box_y end
       if (dz > half_box_z) then dz = dz - box_z end
       if (dx <-half_box_x) then dx = dx + box_x end
       if (dy <-half_box_y) then dy = dy + box_y end
       if (dz <-half_box_z) then dz = dz + box_z end
       var cutoff2 = parts1[part1].core_part_space.cutoff--TODO:regentlib.max(parts1[part1].core_part_space.cutoff, parts2[part2].core_part_space.cutoff)
       cutoff2 = cutoff2 * cutoff2
       var r2 = dx*dx + dy*dy + dz*dz
       if(r2 <= cutoff2) then
         [kernel_name(rexpr parts1[part1] end, rexpr parts2[part2] end, rexpr r2 end)]
       end
     end
   end
end
return pair_task
end

function create_variable_cutoff_self_task( kernel_name )

local task self_task(parts1 : region(ispace(int1d), part), space : region(ispace(int1d), space_config)) where 
      reads(parts1, space), writes(parts1) do
   var box_x = space[0].dim_x
   var box_y = space[0].dim_y
   var box_z = space[0].dim_z
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   var half_box_z = 0.5 * box_z
   for part1 in parts1.ispace do
     for part2 in parts1.ispace do
       --Compute particle distance
       if(part1 < part2) then
         var dx = parts1[part1].core_part_space.pos_x - parts1[part2].core_part_space.pos_x
         var dy = parts1[part1].core_part_space.pos_y - parts1[part2].core_part_space.pos_y
         var dz = parts1[part1].core_part_space.pos_z - parts1[part2].core_part_space.pos_z
         if (dx > half_box_x) then dx = dx - box_x end
         if (dy > half_box_y) then dy = dy - box_y end
         if (dz > half_box_z) then dz = dz - box_z end
         if (dx <-half_box_x) then dx = dx + box_x end
         if (dy <-half_box_y) then dy = dy + box_y end
         if (dz <-half_box_z) then dz = dz + box_z end
         var cutoff2 = parts1[part1].core_part_space.cutoff
         cutoff2 = cutoff2 * cutoff2
         var r2 = dx*dx + dy*dy + dz*dz
         if(r2 <= cutoff2) then
           [kernel_name(rexpr parts1[part1] end, rexpr parts1[part2] end, rexpr r2 end)]
         end
       end
     end 
  end
end

return self_task
end



function create_symmetric_pairtask_runner( kernel_name )

local cell_pair_task = generate_symmetric_pairwise_task( kernel_name )
local cell_self_task = generate_symmetric_self_task( kernel_name )


local task run_symmetric_pairwise_task_code( particles: region(ispace(int1d), part), cell_space : partition(disjoint, particles , ispace(int3d)), space : region(ispace(int1d), space_config))
    where reads(particles, space), writes(particles) do

    --Do all cell2s in the positive direction
    --Not optimised, it does all cell pairs in the domain, doesn't check
    --cutoff radii or anything.
    for cell1 in cell_space.colors do
        cell_self_task(cell_space[cell1], space)
        for cell2 in cell_space.colors do
          if(cell2.x > cell1.x ) then
            cell_pair_task( cell_space[cell1], cell_space[cell2], space )
          elseif( cell2.x == cell1.x and cell2.y > cell1.y ) then
            cell_pair_task( cell_space[cell1], cell_space[cell2], space )
          elseif( cell2.x == cell1.x and cell2.y == cell1.y and cell2.z > cell1.z) then
            cell_pair_task( cell_space[cell1], cell_space[cell2], space )
          end
        end
    end
end
return run_symmetric_pairwise_task_code
end

function create_asymmetric_pairwise_runner( kernel_name )

local cell_pair_task = generate_asymmetric_pairwise_task( kernel_name )
local cell_self_task = generate_asymmetric_self_task( kernel_name )

local task run_asymmetric_pairwise_task_code( particles: region(ispace(int1d), part), cell_space : partition(disjoint, particles , ispace(int3d)), space : region(ispace(int1d), space_config))
    where reads(particles, space), writes(particles) do

    --Not optimised, it does all cell pairs in the domain, doesn't check
    --cutoff radii or anything.
    for cell1 in cell_space.colors do
      cell_self_task(cell_space[cell1], space)
      for cell2 in cell_space.colors do
        if(cell1 ~= cell2) then
          cell_pair_task(cell_space[cell1], cell_space[cell2], space)
        end
      end
    end


end
return run_asymmetric_pairwise_task_code 
end

-----------------------------------------
--End of Launcher Tasks -----------------
-----------------------------------------



--Generate a task to be executed on every particle in the system
function generate_per_part_task( kernel_name )

local task pairwise_task(parts1 : region(ispace(int1d),part), space : region(ispace(int1d), space_config)) where
   reads(parts1, space), writes(parts1) do
   for part1 in parts1.ispace do
         [kernel_name(rexpr parts1[part1] end, rexpr space[0] end)]
   end
end
return pairwise_task
end

function run_per_particle_task( kernel_name )

local per_part_task = generate_per_part_task( kernel_name )
local task run_per_particle_task_code( particles: region(ispace(int1d), part), cell_space : partition(disjoint, particles , ispace(int3d)), space : region(ispace(int1d), space_config) )
    where reads(particles, space), writes(particles) do

    --For each cell, call the task!
    for cell1 in cell_space.colors do
       per_part_task(cell_space[cell1], space)
    end
end

return run_per_particle_task_code
end


-----------------------------------------
--End of Per Part Tasks------------------
-----------------------------------------

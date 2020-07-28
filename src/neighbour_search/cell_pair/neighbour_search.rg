import "regent"

require("defaults")
format = require("std/format")

task zero_neighbour_part(particle_region : region(ispace(int1d), part)) where writes(particle_region.neighbour_part_space) do
  fill(particle_region.neighbour_part_space.cell_id, int3d({0,0,0}))
end

task update_cell_partitions(particles : region(ispace(int1d), part), x_cells : int1d, y_cells : int1d, z_cells : int1d) where
  reads(particles.neighbour_part_space.cell_id) do
  var space_parameter = ispace(int3d, {x_cells, y_cells, z_cells}, {0,0,0})
  var cell_partition = partition(particles.neighbour_part_space.cell_id, space_parameter)
  return cell_partition

end

--Generate the classic MD-style symmetric update pairwise task.
--This function assumes the cutoff is the same for both particles
function generate_symmetric_pairwise_task( kernel_name )
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local part1 = regentlib.newsymbol("part1")
local part2 = regentlib.newsymbol("part2")
local r2 = regentlib.newsymbol(double, "r2")

--local interaction = kernel_name(part1, part2, r2)

local task pairwise_task(parts1 : region(ispace(int1d),part), parts2 : region(ispace(int1d),part), space : region(ispace(int1d), space_config)) where 
   reads(parts1, parts2, space), writes(parts1, parts2) do
   var half_box_x = 0.5 * space[0].dim_x
   var half_box_y = 0.5 * space[0].dim_y
   var half_box_z = 0.5 * space[0].dim_z
   for [part1] in [parts1] do
     for [part2] in [parts2] do
       --Compute particle distance
       var dx = [parts1][part1].core_part_space.pos_x - [parts2][part2].core_part_space.pos_x
       var dy = [parts1][part1].core_part_space.pos_y - [parts2][part2].core_part_space.pos_y
       var dz = [parts1][part1].core_part_space.pos_z - [parts2][part2].core_part_space.pos_z
       if (dx > half_box_x) then dx = dx - half_box_x end
       if (dy > half_box_y) then dy = dy - half_box_y end
       if (dz > half_box_z) then dz = dz - half_box_z end
       if (dx <-half_box_x) then dx = dx + half_box_x end
       if (dy <-half_box_y) then dy = dy + half_box_y end
       if (dz <-half_box_z) then dz = dz + half_box_z end
       var cutoff2 = [parts1][part1].core_part_space.cutoff
       cutoff2 = cutoff2 * cutoff2
       var [r2] = dx*dx + dy*dy + dz*dz
       if([r2] <= cutoff2) then
         [kernel_name(rexpr parts1[part1] end, rexpr parts2[part2] end, rexpr r2 end)]
       end
     end
   end
end
return pairwise_task
end

--Functionality added but unclear on use-case right now
function generate_asymmetric_pairwise_task( kernel_name )
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local part1 = regentlib.newsymbol("part1")
local part2 = regentlib.newsymbol("part2")
local r2 = regentlib.newsymbol("r2")

--Asymmetric kernel can only write to part1
--local interaction = kernel_name(part1, part2, r2)

local task pairwise_task(parts1 : region(ispace(int1d),part), parts2 : region(ispace(int1d),part),  space : region(ispace(int1d), space_config)) where 
   reads(parts1, parts2, space), writes(parts1) do
   var half_box_x = 0.5 * space[0].dim_x
   var half_box_y = 0.5 * space[0].dim_y
   var half_box_z = 0.5 * space[0].dim_z
   for [part1] in [parts1] do
     for [part2] in [parts2] do
       --Compute particle distance
       var dx = [parts1][part1].core_part_space.pos_x - [parts2][part2].core_part_space.pos_x
       var dy = [parts1][part1].core_part_space.pos_y - [parts2][part2].core_part_space.pos_y
       var dz = [parts1][part1].core_part_space.pos_z - [parts2][part2].core_part_space.pos_z
       if (dx > half_box_x) then dx = dx - half_box_x end
       if (dy > half_box_y) then dy = dy - half_box_x end
       if (dz > half_box_z) then dz = dz - half_box_x end
       if (dx <-half_box_x) then dx = dx + half_box_x end
       if (dy <-half_box_y) then dy = dy + half_box_y end
       if (dz <-half_box_z) then dz = dz + half_box_z end
       var cutoff2 = [parts1][part1].core_part_space.cutoff
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
function generate_asymmetric_self_task( kernel_name )
--local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
--local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
--local part1 = regentlib.newsymbol("part1")
--local part2 = regentlib.newsymbol("part2")
--local r2 = regentlib.newsymbol(double, "r2")

--local interaction = kernel_name(part1, part2, r2)

local task pairwise_task(parts1 : region(ispace(int1d),part), parts2 : region(ispace(int1d), part), space : region(ispace(int1d), space_config)) where
   reads(parts1, parts2, space), writes(parts1) do
   var half_box_x = 0.5 * space[0].dim_x
   var half_box_y = 0.5 * space[0].dim_y
   var half_box_z = 0.5 * space[0].dim_z
   for part1 in parts1.ispace do
     for part2 in parts2.ispace do
       --Compute particle distance
       var dx = [parts1][part1].core_part_space.pos_x - [parts2][part2].core_part_space.pos_x
       var dy = [parts1][part1].core_part_space.pos_y - [parts2][part2].core_part_space.pos_y
       var dz = [parts1][part1].core_part_space.pos_z - [parts2][part2].core_part_space.pos_z
       if (dx > half_box_x) then dx = dx - half_box_x end
       if (dy > half_box_y) then dy = dy - half_box_y end
       if (dz > half_box_z) then dz = dz - half_box_z end
       if (dx <-half_box_x) then dx = dx + half_box_x end
       if (dy <-half_box_y) then dy = dy + half_box_y end
       if (dz <-half_box_z) then dz = dz + half_box_z end
       var cutoff2 = [parts1][part1].core_part_space.cutoff
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













--Generate a task to be executed on every particle in the system
function generate_per_part_task( kernel_name )
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local part1 = regentlib.newsymbol("part1")

--local update_particle = kernel_name(part1)

local task pairwise_task(parts1 : region(ispace(int1d),part), space : region(ispace(int1d), space_config)) where
   reads(parts1, space), writes(parts1) do
   for [part1] in [parts1] do
         [kernel_name(rexpr parts1[part1] end)]
   end
end
return pairwise_task
end

function run_per_particle_task( task_name )
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local cell1 = regentlib.newsymbol("cell1")
local cell_space = regentlib.newsymbol("cell_space")
local space = regentlib.newsymbol("space")
local task run_per_particle_task_code( particles: region(ispace(int1d), part), cell_space : partition(disjoint, particles , ispace(int3d)), space : region(ispace(int1d), space_config) )
    where reads(particles, space), writes(particles) do

    --For each cell, call the task!
    for [cell1] in [cell_space].colors do
       task_name([cell_space][cell1], [space])
    end
end

return run_per_particle_task_code
end


-----------------------------------------
--End of Per Part Tasks------------------
-----------------------------------------


--Here we have subset tasks. The idea here is that we have a subset of particles in a cell, and want to compute interactions for only those particles
function generate_asymmetric_pairwise_subset_task( kernel_name )
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local part1 = regentlib.newsymbol("part1")
local part2 = regentlib.newsymbol("part2")
local r2 = regentlib.newsymbol(double, "r2")

--local interaction = kernel_name(part1, part2, r2)
--parts1 is a subset
local task pairwise_task(parts1 : region(ispace(int1d),part), parts2 : region(ispace(int1d),part), space : region(ispace(int1d), space_config)) where
   reads(parts1, parts2, space), writes(parts1) do
   var half_box_x = 0.5 * space[0].dim_x
   var half_box_y = 0.5 * space[0].dim_y
   var half_box_z = 0.5 * space[0].dim_z
   for [part1] in [parts1] do
     for [part2] in [parts2] do
       --Compute particle distance
       var dx = [parts1][part1].core_part_space.pos_x - [parts2][part2].core_part_space.pos_x
       var dy = [parts1][part1].core_part_space.pos_y - [parts2][part2].core_part_space.pos_y
       var dz = [parts1][part1].core_part_space.pos_z - [parts2][part2].core_part_space.pos_z
       if (dx > half_box_x) then dx = dx - half_box_x end
       if (dy > half_box_y) then dy = dy - half_box_y end
       if (dz > half_box_z) then dz = dz - half_box_z end
       if (dx <-half_box_x) then dx = dx + half_box_x end
       if (dy <-half_box_y) then dy = dy + half_box_y end
       if (dz <-half_box_z) then dz = dz + half_box_z end
       var cutoff2 = [parts1][part1].core_part_space.cutoff
       cutoff2 = cutoff2 * cutoff2
       var [r2] = dx*dx + dy*dy + dz*dz
       if([r2] <= cutoff2) then
         [kernel_name(rexpr parts1[part1] end, rexpr parts2[part2] end, rexpr r2 end)]
       end
     end
   end
end
return pairwise_task
end

--This is for if the subset and the cell are overlapping (i.e. the subset is a subset of this cell
function generate_asymmetric_pairwise_subset_task_self( kernel_name )

local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local part1 = regentlib.newsymbol("part1")
local part2 = regentlib.newsymbol("part2")
local r2 = regentlib.newsymbol(double, "r2")

--local interaction = kernel_name(part1, part2, r2)
--parts1 is a subset
local task pairwise_task(parts1 : region(ispace(int1d),part), parts2 : region(ispace(int1d),part), space : region(ispace(int1d), space_config)) where
   reads(parts1, parts2, space), writes(parts1) do
   var half_box_x = 0.5 * space[0].dim_x
   var half_box_y = 0.5 * space[0].dim_y
   var half_box_z = 0.5 * space[0].dim_z
   for [part1] in [parts1] do
     for [part2] in [parts2] do
       --Compute particle distance
       if part1.core_part_space.id ~= part2.core_part_space.id then
         var dx = [parts1][part1].core_part_space.pos_x - [parts2][part2].core_part_space.pos_x
         var dy = [parts1][part1].core_part_space.pos_y - [parts2][part2].core_part_space.pos_y
         var dz = [parts1][part1].core_part_space.pos_z - [parts2][part2].core_part_space.pos_z
         if (dx > half_box_x) then dx = dx - half_box_x end
         if (dy > half_box_y) then dy = dy - half_box_y end
         if (dz > half_box_z) then dz = dz - half_box_z end
         if (dx <-half_box_x) then dx = dx + half_box_x end
         if (dy <-half_box_y) then dy = dy + half_box_y end
         if (dz <-half_box_z) then dz = dz + half_box_z end
         var cutoff2 = [parts1][part1].core_part_space.cutoff
         cutoff2 = cutoff2 * cutoff2
         var [r2] = dx*dx + dy*dy + dz*dz
         --TODO: Remove this assert
         regentlib.assert(r2 ~= 0.0, "r2 is 0, presume self-interaction of a particle..." )
         if([r2] <= cutoff2) then
           [kernel_name(rexpr parts1[part1] end, rexpr parts2[part2] end, rexpr r2 end)]
         end
       end
     end
   end
end
return pairwise_task
end

function run_asymmetric_subset_task( kernel_name )
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local cell1 = regentlib.newsymbol("cell1")
local cell2 = regentlib.newsymbol("cell2")
local cell_space = regentlib.newsymbol("cell_space")
local space = regentlib.newsymbol("space")

local self_task = generate_asymmetric_pairwise_subset_task_self( kernel_name )
local pair_task = generate_asymmetric_pairwise_subset_task( kernel_name )

local task run_asym_subset_task(subset: region(ispace(int1d),part), particles: region(ispace(int1d), part), cell_space : partition(disjoint, particles , ispace(int3d)), space : region(ispace(int1d), space_config))
    where reads(subset,particles, space), writes(subset) do
    var local_cell : int3d
    for s in subset do
      local_cell = subset[s].neighbour_part_space.cell_id
    end
    for cell in cell_space.colors do
      if cell == local_cell then
         self_task(subset, cell_space[cell], space)
      else
         pair_task(subset, cell_space[cell], space)
      end

    end

end

return run_asym_subset_task
end


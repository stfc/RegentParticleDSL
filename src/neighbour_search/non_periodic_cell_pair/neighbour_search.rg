import "regent"

require("defaults")

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
local r2 = regentlib.newsymbol("r2")

local interaction = kernel_name(part1, part2, r2)

local task pairwise_task(parts1 : region(ispace(int1d),part), parts2 : region(ispace(int1d),part), space : region(ispace(int1d), space_config)) where 
   reads(parts1, parts2, space), writes(parts1, parts2) do
   for [part1] in [parts1] do
     for [part2] in [parts2] do
       --Compute particle distance
       var dx = [parts1][part1].core_part_space.pos_x - [parts2][part2].core_part_space.pos_x
       var dy = [parts1][part1].core_part_space.pos_y - [parts2][part2].core_part_space.pos_y
       var dz = [parts1][part1].core_part_space.pos_z - [parts2][part2].core_part_space.pos_z
       var cutoff2 = [parts1][part1].core_part_space.cutoff
       cutoff2 = cutoff2 * cutoff2
       var r2 = dx*dx + dy*dy + dz*dz
       if(r2 <= cutoff2) then
         [interaction]
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
local interaction = kernel_name(part1, part2, r2)

local task pairwise_task(parts1 : region(ispace(int1d),part), parts2 : region(ispace(int1d),part),  space : region(ispace(int1d), space_config)) where 
   reads(parts1, parts2, space), writes(parts1) do
   for [part1] in [parts1] do
     for [part2] in [parts2] do
       --Compute particle distance
       var dx = [parts1][part1].pos_x - [parts2][part2].pos_x
       var dy = [parts1][part1].pos_y - [parts2][part2].pos_y
       var dz = [parts1][part1].pos_z - [parts2][part2].pox_z
       var cutoff2 = [parts1][part1].core_part_space.cutoff
       cutoff2 = cutoff2 * cutoff2
       var r2 = dx*dx + dy*dy + dz*dz
       if(r2 <= cutoff2) then
         [interaction]
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

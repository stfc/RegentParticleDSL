-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

require("defaults")
require("src/neighbour_search/cell_pair_v2/cell")
compute_privileges = require("src/utils/compute_privilege")
format = require("std/format")
string_to_field_path = require("src/utils/string_to_fieldpath")

task zero_neighbour_part(particle_region : region(ispace(int1d), part)) where writes(particle_region.neighbour_part_space) do
  fill(particle_region.neighbour_part_space.cell_id, int3d({0,0,0}))
end

__demand(__inline)
task update_cell_partitions(particles : region(ispace(int1d), part), config : region(ispace(int1d), config_type)) where
  reads(particles.neighbour_part_space,particles.core_part_space, config),
  writes(config.neighbour_config, particles.neighbour_part_space) do
  initialise_cells(config, particles)
  particles_to_cell_launcher(particles, config)
  var space_parameter = ispace(int3d, {config[0].neighbour_config.x_cells, config[0].neighbour_config.y_cells, config[0].neighbour_config.z_cells}, {0,0,0})
  var cell_partition = partition(particles.neighbour_part_space.cell_id, space_parameter)
  return cell_partition

end

local abs = regentlib.fabs(double)
local ceil = regentlib.ceil(double)

--Generate the classic MD-style symmetric update pairwise task.
--This function assumes the cutoff is the same for both particles
function generate_symmetric_pairwise_task( kernel_name, read1, read2, write1, write2 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local read1_privs = terralib.newlist()
local read2_privs = terralib.newlist()
local write1_privs = terralib.newlist()
local write2_privs = terralib.newlist()

for _, v in pairs(read1) do
  read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v))) 
end
for _, v in pairs(read2) do
  read2_privs:insert( regentlib.privilege(regentlib.reads, parts2, string_to_field_path.get_field_path(v))) 
end
for _, v in pairs(write1) do
  write1_privs:insert( regentlib.privilege(regentlib.writes, parts1, string_to_field_path.get_field_path(v))) 
end
for _, v in pairs(write2) do
  write2_privs:insert( regentlib.privilege(regentlib.writes, parts2, string_to_field_path.get_field_path(v))) 
end

local task pairwise_task([parts1], [parts2], config : region(ispace(int1d), config_type))  
  where [read1_privs], [read2_privs], [write1_privs], [write2_privs], reads(config), reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}),
  reads(parts2.core_part_space.{pos_x, pos_y, pos_z, cutoff}) do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var box_z = config[0].space.dim_z
   var half_box_x = 0.5 * box_x 
   var half_box_y = 0.5 * box_y
   var half_box_z = 0.5 * box_z
   for part1 in [parts1].ispace do
     for part2 in [parts2].ispace do
       --Compute particle distance
       var dx = [parts1][part1].core_part_space.pos_x - [parts2][part2].core_part_space.pos_x
       var dy = [parts1][part1].core_part_space.pos_y - [parts2][part2].core_part_space.pos_y
       var dz = [parts1][part1].core_part_space.pos_z - [parts2][part2].core_part_space.pos_z
       if (dx > half_box_x) then dx = dx - box_x end
       if (dy > half_box_y) then dy = dy - box_y end
       if (dz > half_box_z) then dz = dz - box_z end
       if (dx <-half_box_x) then dx = dx + box_x end
       if (dy <-half_box_y) then dy = dy + box_y end
       if (dz <-half_box_z) then dz = dz + box_z end
       var cutoff2 = regentlib.fmax([parts1][part1].core_part_space.cutoff, [parts2][part2].core_part_space.cutoff)
       cutoff2 = cutoff2 * cutoff2
       var r2 = dx*dx + dy*dy + dz*dz
       if(r2 <= cutoff2) then
         [kernel_name(rexpr [parts1][part1] end, rexpr [parts2][part2] end, rexpr r2 end)]
       end
     end
   end
end
return pairwise_task
end

--Functionality added but unclear on use-case right now
function generate_asymmetric_pairwise_task( kernel_name, read1, read2, write1 )
--Asymmetric kernel can only write to part1
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local read1_privs = terralib.newlist()
local read2_privs = terralib.newlist()
local write1_privs = terralib.newlist()

for _, v in pairs(read1) do
  read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
end
for _, v in pairs(read2) do
  read2_privs:insert( regentlib.privilege(regentlib.reads, parts2, string_to_field_path.get_field_path(v)))
end
for _, v in pairs(write1) do
  write1_privs:insert( regentlib.privilege(regentlib.writes, parts1, string_to_field_path.get_field_path(v)))
end

local task pairwise_task([parts1], [parts2],  config : region(ispace(int1d), config_type))
  where [read1_privs], [read2_privs], [write1_privs], reads(config), reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}), 
  reads(parts2.core_part_space.{pos_x, pos_y, pos_z, cutoff}) do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var box_z = config[0].space.dim_z
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   var half_box_z = 0.5 * box_z
   for part1 in [parts1].ispace do
     for part2 in [parts2].ispace do
       --Compute particle distance
       var dx = [parts1][part1].core_part_space.pos_x - [parts2][part2].core_part_space.pos_x
       var dy = [parts1][part1].core_part_space.pos_y - [parts2][part2].core_part_space.pos_y
       var dz = [parts1][part1].core_part_space.pos_z - [parts2][part2].core_part_space.pos_z
       if (dx > half_box_x) then dx = dx - box_x end
       if (dy > half_box_y) then dy = dy - box_y end
       if (dz > half_box_z) then dz = dz - box_z end
       if (dx <-half_box_x) then dx = dx + box_x end
       if (dy <-half_box_y) then dy = dy + box_y end
       if (dz <-half_box_z) then dz = dz + box_z end
       var cutoff2 = [parts1][part1].core_part_space.cutoff
       cutoff2 = cutoff2 * cutoff2
       var r2 = dx*dx + dy*dy + dz*dz
       if(r2 <= cutoff2) then
         [kernel_name(rexpr [parts1][part1] end, rexpr [parts2][part2] end, rexpr r2 end)]
       end
     end
   end
end
return pairwise_task
end

-----------------------------------------
--End of Pair Tasks ---------------------
-----------------------------------------


--Generate a self task
--This function assumes the cutoff is the same for both particles
function generate_symmetric_self_task( kernel_name, read1, read2, write1, write2 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local read1_privs = terralib.newlist()
local write1_privs = terralib.newlist()

for _, v in pairs(read1) do
  read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
end
for _, v in pairs(read2) do
  read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
end
for _, v in pairs(write1) do
  write1_privs:insert( regentlib.privilege(regentlib.writes, parts1, string_to_field_path.get_field_path(v)))
end
for _, v in pairs(write2) do
  write1_privs:insert( regentlib.privilege(regentlib.writes, parts1, string_to_field_path.get_field_path(v)))
end

local task self_task([parts1], config : region(ispace(int1d),config_type)) where
  [read1_privs], [write1_privs], reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}), reads(config) do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var box_z = config[0].space.dim_z
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   var half_box_z = 0.5 * box_z
   for part1 in [parts1].ispace do
     for part2 in [parts1].ispace do
       --Compute particle distance
       if(part1 < part2) then
         var dx = [parts1][part1].core_part_space.pos_x - [parts1][part2].core_part_space.pos_x
         var dy = [parts1][part1].core_part_space.pos_y - [parts1][part2].core_part_space.pos_y
         var dz = [parts1][part1].core_part_space.pos_z - [parts1][part2].core_part_space.pos_z
         if (dx > half_box_x) then dx = dx - box_x end
         if (dy > half_box_y) then dy = dy - box_y end
         if (dz > half_box_z) then dz = dz - box_z end
         if (dx <-half_box_x) then dx = dx + box_x end
         if (dy <-half_box_y) then dy = dy + box_y end
         if (dz <-half_box_z) then dz = dz + box_z end
         var cutoff2 = regentlib.fmax([parts1][part1].core_part_space.cutoff, [parts1][part2].core_part_space.cutoff)
         cutoff2 = cutoff2 * cutoff2
         var r2 = dx*dx + dy*dy + dz*dz
         if(r2 <= cutoff2) then
           [kernel_name(rexpr [parts1][part1] end, rexpr [parts1][part2] end, rexpr r2 end)]
         end
       end
     end
   end

end
return self_task
end

--Generate a self task
--This function assumes the cutoff of only the updated part is relevant
function generate_asymmetric_self_task( kernel_name, read1, read2, write1 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local read1_privs = terralib.newlist()
local write1_privs = terralib.newlist()

for _, v in pairs(read1) do
  read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
end
for _, v in pairs(read2) do
  read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
end
for _, v in pairs(write1) do
  write1_privs:insert( regentlib.privilege(regentlib.writes, parts1, string_to_field_path.get_field_path(v)))
end

local task self_task([parts1], config : region(ispace(int1d), config_type)) where
   [read1_privs], [write1_privs], reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}), reads(config) do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var box_z = config[0].space.dim_z
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   var half_box_z = 0.5 * box_z
   for part1 in [parts1].ispace do
     for part2 in [parts1].ispace do
       --Compute particle distance
       if(part1 ~= part2) then
         var dx = [parts1][part1].core_part_space.pos_x - [parts1][part2].core_part_space.pos_x
         var dy = [parts1][part1].core_part_space.pos_y - [parts1][part2].core_part_space.pos_y
         var dz = [parts1][part1].core_part_space.pos_z - [parts1][part2].core_part_space.pos_z
         if (dx > half_box_x) then dx = dx - box_x end
         if (dy > half_box_y) then dy = dy - box_y end
         if (dz > half_box_z) then dz = dz - box_z end
         if (dx <-half_box_x) then dx = dx + box_x end
         if (dy <-half_box_y) then dy = dy + box_y end
         if (dz <-half_box_z) then dz = dz + box_z end
         var cutoff2 = [parts1][part1].core_part_space.cutoff
         cutoff2 = cutoff2 * cutoff2
         var r2 = dx*dx + dy*dy + dz*dz
         if(r2 <= cutoff2) then
           [kernel_name(rexpr [parts1][part1] end, rexpr [parts1][part2] end, rexpr r2 end)]
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

local task pair_task(parts1 : region(ispace(int1d), part), parts2 : region(ispace(int1d), part), config : region(ispace(int1d), config_type)) where
     reads(parts1, parts2, config), writes( parts1, parts2) do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var box_z = config[0].space.dim_z
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

local task self_task(parts1 : region(ispace(int1d), part), config : region(ispace(int1d), config_type)) where 
      reads(parts1, config), writes(parts1) do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var box_z = config[0].space.dim_z
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
         var cutoff2 = parts1[part1].core_part_space.cutoff--TODO:regentlib.max(parts1[part1].core_part_space.cutoff, parts2[part2].core_part_space.cutoff)
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

__demand(__inline)
task cell_greater_equal(cell1 : int3d, cell2 : int3d) : bool
  var returnval : bool = false
  if (cell1.x > cell2.x or (cell1.x == cell2.x and cell1.y > cell2.y) or( cell1.x == cell2.x and cell1.y == cell2.y and cell1.z > cell2.z)) then
    returnval =  true
  end
  return returnval
end

function create_symmetric_pairwise_runner( kernel_name, config, cell_space )
local read1, read2, write1, write2 = compute_privileges.two_region_privileges( kernel_name )
local cell_pair_task = generate_symmetric_pairwise_task( kernel_name, read1, read2, write1, write2 )
local cell_self_task = generate_symmetric_self_task( kernel_name, read1, read2, write1, write2 )

local symmetric = rquote
    --Do all cell2s in the positive direction
    var cutoff2 = config[0].neighbour_config.max_cutoff * config[0].neighbour_config.max_cutoff
    var x_count = config[0].neighbour_config.x_cells
    var y_count = config[0].neighbour_config.y_cells
    var z_count = config[0].neighbour_config.z_cells

    --Compute cell radii
    var cutoff = config[0].neighbour_config.max_cutoff
    var x_radii : int = ceil( cutoff / config[0].neighbour_config.cell_dim_x )
    var y_radii : int = ceil( cutoff / config[0].neighbour_config.cell_dim_y )
    var z_radii : int = ceil( cutoff / config[0].neighbour_config.cell_dim_z )
    for cell1 in cell_space.colors do
        cell_self_task(cell_space[cell1], config)
        --Loops non inclusive, positive only direction.
        for x = -x_radii, x_radii+1 do
          for y = -y_radii, y_radii+1 do
            for z = -z_radii, z_radii + 1 do
              if(not (x == 0 and y == 0 and z == 0) ) then
                var cell2_x = cell1.x + x
                var cell2_y = cell1.y + y
                var cell2_z = cell1.z + z
                if cell2_x < 0 then
                  cell2_x = cell2_x + x_count
                end
                if cell2_y < 0 then
                  cell2_y = cell2_y + y_count
                end
                if cell2_z < 0 then
                  cell2_z = cell2_z + z_count
                end
                var cell2 : int3d = int3d({ (cell2_x)%x_count, (cell2_y)%y_count, (cell2_z)%z_count })
                --Weird if statement to handle max_cutoff >= half the boxsize
                if( cell_greater_equal(cell1, cell2) 
                or ((x < 0 and cell2.x <= cell1.x + x_radii and cell2.x > cell1.x) or
                                                        (y < 0 and cell2.y <= cell1.y + y_radii and cell2.y > cell1.y) or
                                                        (z < 0 and cell2.z <= cell1.z + z_radii and cell2.z > cell1.z))
                or ((x > 0 and cell2.x >= cell1.x - x_radii and cell2.x < cell1.x) or (y > 0 and cell2.y >= cell1.y - y_radii and cell2.y < cell1.y)
                    or (z > 0 and cell2.z >= cell1.z - z_radii and cell2.z < cell1.z) )
                ) then
                else
                  --symmetric
                  cell_pair_task(cell_space[cell1], cell_space[cell2], config)
                end
              end
            end
          end
        end
    end
end
return symmetric

end

function create_asymmetric_pairwise_runner( kernel_name, config, cell_space )
local read1, read2, write1, write2 = compute_privileges.two_region_privileges( kernel_name )
--While write2 is computed, asymmetric kernels are not allowed to write to write2
local cell_pair_task = generate_asymmetric_pairwise_task( kernel_name, read1, read2, write1 )
local cell_self_task = generate_asymmetric_self_task( kernel_name, read1, read2, write1 )

local asymmetric = rquote
    --Do all cell2s in the positive direction
    var cutoff2 = config[0].neighbour_config.max_cutoff * config[0].neighbour_config.max_cutoff
    var x_count = config[0].neighbour_config.x_cells
    var y_count = config[0].neighbour_config.y_cells
    var z_count = config[0].neighbour_config.z_cells

    --Compute cell radii
    var cutoff = config[0].neighbour_config.max_cutoff
    var x_radii : int = ceil( cutoff / config[0].neighbour_config.cell_dim_x )
    var y_radii : int = ceil( cutoff / config[0].neighbour_config.cell_dim_y )
    var z_radii : int = ceil( cutoff / config[0].neighbour_config.cell_dim_z )
    for cell1 in cell_space.colors do
        cell_self_task(cell_space[cell1], config)
        --Loops non inclusive, positive only direction.
        for x = -x_radii, x_radii+1 do
          for y = -y_radii, y_radii+1 do
            for z = -z_radii, z_radii + 1 do
              if(not (x == 0 and y == 0 and z == 0) ) then
                var cell2_x = cell1.x + x
                var cell2_y = cell1.y + y
                var cell2_z = cell1.z + z
                if cell2_x < 0 then
                  cell2_x = cell2_x + x_count
                end
                if cell2_y < 0 then
                  cell2_y = cell2_y + y_count
                end
                if cell2_z < 0 then
                  cell2_z = cell2_z + z_count
                end
                var cell2 : int3d = int3d({ (cell2_x)%x_count, (cell2_y)%y_count, (cell2_z)%z_count })
                --if statement to handle max_cutoff >= half the boxsize
                    --Handle radius overlap
                if( cell_greater_equal(cell1, cell2) 
                or ((x < 0 and cell2.x <= cell1.x + x_radii and cell2.x > cell1.x) or
                                                        (y < 0 and cell2.y <= cell1.y + y_radii and cell2.y > cell1.y) or
                                                        (z < 0 and cell2.z <= cell1.z + z_radii and cell2.z > cell1.z))
                or ((x > 0 and cell2.x >= cell1.x - x_radii and cell2.x < cell1.x) or (y > 0 and cell2.y >= cell1.y - y_radii and cell2.y < cell1.y)
                    or (z > 0 and cell2.z >= cell1.z - z_radii and cell2.z < cell1.z) )
                ) then
                 else
                  --Asymmetric, launch tasks both ways
                  cell_pair_task(cell_space[cell1], cell_space[cell2], config)
                  cell_pair_task(cell_space[cell2], cell_space[cell1], config)
                end
              end
            end
          end
        end
    end

end
return asymmetric
end


-----------------------------------------
--End of Launcher Tasks -----------------
-----------------------------------------



--Generate a task to be executed on every particle in the system
function generate_per_part_task( kernel_name, read1, write1 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local read1_privs = terralib.newlist()
local write1_privs = terralib.newlist()

for _, v in pairs(read1) do
  read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
end
for _, v in pairs(write1) do
  write1_privs:insert( regentlib.privilege(regentlib.writes, parts1, string_to_field_path.get_field_path(v)))
end
local task pairwise_task([parts1], config : region(ispace(int1d), config_type)) where
   [read1_privs], [write1_privs], reads(config) do
   for part1 in [parts1].ispace do
         [kernel_name(rexpr [parts1][part1] end, rexpr config[0] end)]
   end
end
return pairwise_task
end

function generate_per_part_task_bool_return ( kernel_name )

local task per_part_bool_task(parts1 : region(ispace(int1d),part), config : region(ispace(int1d), config_type)) : bool where
   reads(parts1, config), writes(parts1) do
   var return_val : bool = false
   for part1 in parts1.ispace do
         [kernel_name(rexpr parts1[part1] end, rexpr config[0] end, rexpr return_val end)]
   end
   return return_val
end
return per_part_bool_task
end


function run_per_particle_task( kernel_name, config, cell_space )

local read1, read2, write1, write2 = compute_privileges.two_region_privileges(kernel_name)
local per_part_task = generate_per_part_task( kernel_name, read1, write1 )
local runner = rquote

    --For each cell, call the task!
    for cell1 in cell_space.colors do
       per_part_task(cell_space[cell1], config)
    end
end

return runner
end

-----------------------------------------
--End of Per Part Tasks------------------
-----------------------------------------

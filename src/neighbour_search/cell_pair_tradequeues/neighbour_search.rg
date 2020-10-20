-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

require("defaults")
require("src/neighbour_search/cell_pair_tradequeues/cell")
local compute_privileges = require("src/utils/compute_privilege")
local format = require("std/format")
local string_to_field_path = require("src/utils/string_to_fieldpath")


local abs = regentlib.fabs(double)
local ceil = regentlib.ceil(double)

--This function assumes the cutoff is the same for both particles
function generate_symmetric_pairwise_task( kernel_name, read1, read2, write1, write2 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local read1_privs = terralib.newlist()
local read2_privs = terralib.newlist()
local write1_privs = terralib.newlist()
local write2_privs = terralib.newlist()
local update_neighbours = false

for _, v in pairs(read1) do
  read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
end
for _, v in pairs(read2) do
  read2_privs:insert( regentlib.privilege(regentlib.reads, parts2, string_to_field_path.get_field_path(v)))
end
for _, v in pairs(write1) do
  write1_privs:insert( regentlib.privilege(regentlib.writes, parts1, string_to_field_path.get_field_path(v)))
  if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" then
    update_neighbours = true
  end
end
for _, v in pairs(write2) do
  write2_privs:insert( regentlib.privilege(regentlib.writes, parts2, string_to_field_path.get_field_path(v)))
  if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" then
    update_neighbours = true
  end
end

local task pairwise_task([parts1], [parts2], config : region(ispace(int1d), config_type))
  where [read1_privs], [read2_privs], [write1_privs], [write2_privs], reads(config), reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}),
  reads(parts2.core_part_space.{pos_x, pos_y, pos_z, cutoff}), reads(parts1.neighbour_part_space._valid), reads(parts2.neighbour_part_space._valid) do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var box_z = config[0].space.dim_z
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   var half_box_z = 0.5 * box_z
   --For the tradequeue implementation, we only care about the validity of the particles
   for part1 in [parts1].ispace do
     if [parts1][part1].neighbour_part_space._valid then
       for part2 in [parts2].ispace do
         if [parts2][part2].neighbour_part_space._valid then
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
   end
end
return pairwise_task, update_neighbours
end


function generate_asymmetric_pairwise_task( kernel_name, read1, read2, write1 )
--Asymmetric kernel can only write to part1
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local read1_privs = terralib.newlist()
local read2_privs = terralib.newlist()
local write1_privs = terralib.newlist()
local update_neighbours = false

for _, v in pairs(read1) do
  read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
end
for _, v in pairs(read2) do
  read2_privs:insert( regentlib.privilege(regentlib.reads, parts2, string_to_field_path.get_field_path(v)))
end
for _, v in pairs(write1) do
  write1_privs:insert( regentlib.privilege(regentlib.writes, parts1, string_to_field_path.get_field_path(v)))
  if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" then
    update_neighbours = true
  end
end
local task pairwise_task([parts1], [parts2],  config : region(ispace(int1d), config_type))
  where [read1_privs], [read2_privs], [write1_privs], reads(config), reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}),
  reads(parts2.core_part_space.{pos_x, pos_y, pos_z, cutoff}), reads(parts1.neighbour_part_space._valid), reads(parts2.neighbour_part_space._valid) do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var box_z = config[0].space.dim_z
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   var half_box_z = 0.5 * box_z
   for part1 in [parts1].ispace do
     if [parts1][part1].neighbour_part_space._valid then
       for part2 in [parts2].ispace do
         if [parts2][part2].neighbour_part_space._valid then
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
   end
end
return pairwise_task, update_neighbours
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
local update_neighbours = false

for _, v in pairs(read1) do
  read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
end
for _, v in pairs(read2) do
  read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
end
for _, v in pairs(write1) do
  write1_privs:insert( regentlib.privilege(regentlib.writes, parts1, string_to_field_path.get_field_path(v)))
  if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" then
    update_neighbours = true
  end
end
for _, v in pairs(write2) do
  write1_privs:insert( regentlib.privilege(regentlib.writes, parts1, string_to_field_path.get_field_path(v)))
  if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" then
    update_neighbours = true
  end
end

local task self_task([parts1], config : region(ispace(int1d),config_type)) where
  [read1_privs], [write1_privs], reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}),
                                 reads(parts1.neighbour_part_space._valid), reads(config) do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var box_z = config[0].space.dim_z
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   var half_box_z = 0.5 * box_z
   for part1 in [parts1].ispace do
     if [parts1][part1].neighbour_part_space._valid then
       for part2 in [parts1].ispace do
         if [parts1][part2].neighbour_part_space._valid then
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
   end

end
return self_task, update_neighbours
end


--Generate a self task
--This function assumes the cutoff of only the updated part is relevant
function generate_asymmetric_self_task( kernel_name, read1, read2, write1 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local read1_privs = terralib.newlist()
local write1_privs = terralib.newlist()
local update_neighbours = false

for _, v in pairs(read1) do
  read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
end
for _, v in pairs(read2) do
  read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
end
for _, v in pairs(write1) do
  write1_privs:insert( regentlib.privilege(regentlib.writes, parts1, string_to_field_path.get_field_path(v)))
  if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" then
    update_neighbours = true
  end
end

local task self_task([parts1], config : region(ispace(int1d), config_type)) where
   [read1_privs], [write1_privs], reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}), 
                                  reads(parts1.neighbour_part_space._valid), reads(config) do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var box_z = config[0].space.dim_z
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   var half_box_z = 0.5 * box_z
   for part1 in [parts1].ispace do
     if [parts1][part1].neighbour_part_space._valid then
       for part2 in [parts1].ispace do
         if [parts1][part2].neighbour_part_space._valid then
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
   end
end
return self_task, update_neighbours
end

-----------------------------------------
--End of Self Tasks ---------------------
-----------------------------------------


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
local cell_pair_task, update_neighbours1 = generate_symmetric_pairwise_task( kernel_name, read1, read2, write1, write2 )
local cell_self_task, update_neighbours2 = generate_symmetric_self_task( kernel_name, read1, read2, write1, write2 )
local update_neighbours = update_neighbours1 or update_neighbours2
local update_neighbours_quote = rquote

end
if update_neighbours then
  local temp_variables = {}
  temp_variables.config = config
  update_neighbours_quote = neighbour_init.update_cells(temp_variables)
end

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
    [update_neighbours_quote];
end
return symmetric

end

function create_asymmetric_pairwise_runner( kernel_name, config, cell_space )
local read1, read2, write1, write2 = compute_privileges.two_region_privileges( kernel_name )
--While write2 is computed, asymmetric kernels are not allowed to write to write2
local cell_pair_task, update_neighbours1 = generate_asymmetric_pairwise_task( kernel_name, read1, read2, write1 )
local cell_self_task, update_neighbours2 = generate_asymmetric_self_task( kernel_name, read1, read2, write1 )
local update_neighbours = update_neighbours1 or update_neighbours2
local update_neighbours_quote = rquote

end
if update_neighbours then
  local temp_variables = {}
  temp_variables.config = config
  update_neighbours_quote = neighbour_init.update_cells(temp_variables)
end


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
    [update_neighbours_quote];
end
return asymmetric
end


-----------------------------------------
--End of Launcher Code  -----------------
-----------------------------------------



--Generate a task to be executed on every particle in the system
function generate_per_part_task( kernel_name, read1, write1 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local read1_privs = terralib.newlist()
local write1_privs = terralib.newlist()
local update_neighbours = false
for _, v in pairs(read1) do
  read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
end
--Add the extra privilege required by the search algorithm
read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path("neighbour_part_space._valid")  ))
for _, v in pairs(write1) do
  write1_privs:insert( regentlib.privilege(regentlib.writes, parts1, string_to_field_path.get_field_path(v)))
  if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" then
    update_neighbours = true
  end
end
local task pairwise_task([parts1], config : region(ispace(int1d), config_type)) where
   [read1_privs], [write1_privs], reads(config) do
   for part1 in [parts1].ispace do
     if [parts1][part1].neighbour_part_space._valid then
         [kernel_name(rexpr [parts1][part1] end, rexpr config[0] end)]
     end
   end
end
return pairwise_task, update_neighbours
end

function generate_per_part_task_bool_return ( kernel_name )

local task per_part_bool_task(parts1 : region(ispace(int1d),part), config : region(ispace(int1d), config_type)) : bool where
   reads(parts1, config), writes(parts1) do
   var return_val : bool = false
   for part1 in parts1.ispace do
     if parts1[part1].neighbour_part_space._valid then 
         [kernel_name(rexpr parts1[part1] end, rexpr config[0] end, rexpr return_val end)]
     end
   end
   return return_val
end
return per_part_bool_task
end

function run_per_particle_task( kernel_name, config, cell_space )

local read1, read2, write1, write2 = compute_privileges.two_region_privileges(kernel_name)
local need_tradequeues = false
local per_part_task, update_neighbours = generate_per_part_task( kernel_name, read1, write1 )
local update_neighbours_quote = rquote

end
if update_neighbours then
  local temp_variables = {}
  temp_variables.config = config
  update_neighbours_quote = neighbour_init.update_cells(temp_variables)
end
local runner = rquote

    --For each cell, call the task!
    for cell1 in cell_space.colors do
       per_part_task(cell_space[cell1], config)
    end
   [update_neighbours_quote];
end

return runner
end

-----------------------------------------
--End of Per Part Tasks------------------
-----------------------------------------

-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

require("src/neighbour_search/cell_pair_tradequeues_nonperiod/cell")
local compute_privileges = require("src/utils/compute_privilege")
local format = require("std/format")
local string_to_field_path = require("src/utils/string_to_fieldpath")
local c = regentlib.c

local abs = regentlib.fabs(double)
local ceil = regentlib.ceil(double)


--This function assumes the cutoff is the same for both particles
function generate_symmetric_pairwise_task_multikernel( kernel_list, read1, read2, write1, write2 )
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

local coherences = terralib.newlist()
if update_neighbours then
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts2 ) )
else
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts2 ) )
end

local __demand(__leaf) task pairwise_task([parts1], [parts2], config : region(ispace(int1d), config_type))
  where [read1_privs], [read2_privs], [write1_privs], [write2_privs], reads(config), reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}),
  reads(parts2.core_part_space.{pos_x, pos_y, pos_z, cutoff}), reads(parts1.neighbour_part_space._valid), reads(parts2.neighbour_part_space._valid),
 [coherences] do
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
           var cutoff2 = regentlib.fmax([parts1][part1].core_part_space.cutoff, [parts2][part2].core_part_space.cutoff)
           cutoff2 = cutoff2 * cutoff2
           var r2 = dx*dx + dy*dy + dz*dz
           if(r2 <= cutoff2) then
             [kernel_list:map( function(kernel)
               return kernel(rexpr [parts1][part1] end, rexpr [parts2][part2] end, rexpr r2 end)
             end)];
           end
         end
       end
     end
   end
end
return pairwise_task, update_neighbours
end


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

local coherences = terralib.newlist()
if update_neighbours then
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts2 ) )
else
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts2 ) ) 
end

local __demand(__leaf) task pairwise_task([parts1], [parts2], config : region(ispace(int1d), config_type))
  where [read1_privs], [read2_privs], [write1_privs], [write2_privs], reads(config), reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}),
  reads(parts2.core_part_space.{pos_x, pos_y, pos_z, cutoff}), reads(parts1.neighbour_part_space._valid), reads(parts2.neighbour_part_space._valid),
 [coherences] do
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

function generate_asymmetric_pairwise_task_multikernel( kernel_list, read1, read2, write1 )
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

local coherences = terralib.newlist()
if update_neighbours then
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts2 ) )
else
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts2 ) )
end

local __demand(__leaf) task pairwise_task([parts1], [parts2],  config : region(ispace(int1d), config_type))
  where [read1_privs], [read2_privs], [write1_privs], reads(config), reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}),
  reads(parts2.core_part_space.{pos_x, pos_y, pos_z, cutoff}), reads(parts1.neighbour_part_space._valid), reads(parts2.neighbour_part_space._valid),
  [coherences] do
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
           var cutoff2 = [parts1][part1].core_part_space.cutoff
           cutoff2 = cutoff2 * cutoff2
           var r2 = dx*dx + dy*dy + dz*dz
           if(r2 <= cutoff2) then
             [kernel_list:map( function(kernel)
               return kernel(rexpr [parts1][part1] end, rexpr [parts2][part2] end, rexpr r2 end)
             end)];
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

local coherences = terralib.newlist()
if update_neighbours then
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts2 ) )
else
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts2 ) )
end

local __demand(__leaf) task pairwise_task([parts1], [parts2],  config : region(ispace(int1d), config_type))
  where [read1_privs], [read2_privs], [write1_privs], reads(config), reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}),
  reads(parts2.core_part_space.{pos_x, pos_y, pos_z, cutoff}), reads(parts1.neighbour_part_space._valid), reads(parts2.neighbour_part_space._valid),
  [coherences] do
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

function generate_symmetric_self_task_multikernel( kernel_list, read1, read2, write1, write2 )
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

local coherences = terralib.newlist()
if update_neighbours then
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
else
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
end

local __demand(__leaf) task self_task([parts1], config : region(ispace(int1d),config_type)) where
  [read1_privs], [write1_privs], reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}),
                                 reads(parts1.neighbour_part_space._valid), reads(config),
   [coherences] do
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
             var cutoff2 = regentlib.fmax([parts1][part1].core_part_space.cutoff, [parts1][part2].core_part_space.cutoff)
             cutoff2 = cutoff2 * cutoff2
             var r2 = dx*dx + dy*dy + dz*dz
             if(r2 <= cutoff2) then
               [kernel_list:map( function(kernel)
                 return kernel(rexpr [parts1][part1] end, rexpr [parts1][part2] end, rexpr r2 end)
               end)];
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

local coherences = terralib.newlist()
if update_neighbours then
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
else
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
end

local __demand(__leaf) task self_task([parts1], config : region(ispace(int1d),config_type)) where
  [read1_privs], [write1_privs], reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}),
                                 reads(parts1.neighbour_part_space._valid), reads(config),
   [coherences] do
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

function generate_asymmetric_self_task_multikernel( kernel_list, read1, read2, write1 )
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

local coherences = terralib.newlist()
if update_neighbours then
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
else
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
end

local __demand(__leaf) task self_task([parts1], config : region(ispace(int1d), config_type)) where
   [read1_privs], [write1_privs], reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}),
                                  reads(parts1.neighbour_part_space._valid), reads(config),
   [coherences] do
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
             var cutoff2 = [parts1][part1].core_part_space.cutoff
             cutoff2 = cutoff2 * cutoff2
             var r2 = dx*dx + dy*dy + dz*dz
             if(r2 <= cutoff2) then
               [kernel_list:map( function(kernel)
                 return kernel(rexpr [parts1][part1] end, rexpr [parts1][part2] end, rexpr r2 end)
               end)];
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

local coherences = terralib.newlist()
if update_neighbours then
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
else
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
end

local __demand(__leaf) task self_task([parts1], config : region(ispace(int1d), config_type)) where
   [read1_privs], [write1_privs], reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}), 
                                  reads(parts1.neighbour_part_space._valid), reads(config),
   [coherences] do
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
                --No periodic wrapping
                if (cell2_x >= 0 and cell2_x < x_count) 
                   and (cell2_y >= 0 and cell2_y < y_count)
                   and (cell2_z >= 0 and cell2_z < z_count) then
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
                --No periodic wrapping
                if (cell2_x >= 0 and cell2_x < x_count) 
                   and (cell2_y >= 0 and cell2_y < y_count)
                   and (cell2_z >= 0 and cell2_z < z_count) then
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
    [update_neighbours_quote];
end
return asymmetric
end


-----------------------------------------
--End of Launcher Code  -----------------
-----------------------------------------

--Generate a task to be executed on every particle in the system, using a list of kernels
function generate_per_part_task_multikernel( kernel_list, read1, write1 )
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
local coherences = terralib.newlist()
if update_neighbours then
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
else
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
end

local __demand(__leaf) task pairwise_task([parts1], config : region(ispace(int1d), config_type)) where
   [read1_privs], [write1_privs], reads(config), [coherences] do
   for part1 in [parts1].ispace do
     if [parts1][part1].neighbour_part_space._valid then
       [kernel_list:map( function(kernel) 
         return kernel(rexpr [parts1][part1] end, rexpr config[0] end)
       end)];
     end
   end
end
return pairwise_task, update_neighbours

end


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
local coherences = terralib.newlist()
if update_neighbours then
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
else
  coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
end
local __demand(__leaf) task pairwise_task([parts1], config : region(ispace(int1d), config_type)) where
   [read1_privs], [write1_privs], reads(config), [coherences] do
   for part1 in [parts1].ispace do
     if [parts1][part1].neighbour_part_space._valid then
         [kernel_name(rexpr [parts1][part1] end, rexpr config[0] end)]
     end
   end
end
return pairwise_task, update_neighbours
end

function generate_per_part_task_bool_return ( kernel_name )

local __demand(__leaf) task per_part_bool_task(parts1 : region(ispace(int1d),part), config : region(ispace(int1d), config_type)) : bool where
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

-----------------------------------------
--Multikernel Launchers------------------
-----------------------------------------

-- ... should contain a list of symmetric kernel functions
function create_symmetric_pairwise_runner_multikernel( config, cell_space, ... )
  local read1 = terralib.newlist()
  local read2 = terralib.newlist()
  local write1 = terralib.newlist()
  local write2 = terralib.newlist()
  if select("#", ...) == 0 then
    print("Attempting to make a 0 kernel task")
    return nil
  end
  local kernels = terralib.newlist()
  local hash_r1 = {}
  local hash_r2 = {}
  local hash_w1 = {}
  local hash_w2 = {}
  for i=1, select("#", ...) do
    kernels:insert(select(i, ...))
    local temp_r1, temp_r2, temp_w1, temp_w2 = compute_privileges.two_region_privileges(select(i, ...))
    --Merge the read/writes for this kernel with previous ones, keeping uniqueness
    for _,v in pairs(temp_r1) do
      if( not hash_r1[v]) then
        read1:insert(v)
        hash_r1[v] = true
      end
    end
    for _,v in pairs(temp_r2) do
      if( not hash_r2[v]) then
        read2:insert(v)
        hash_r2[v] = true
      end
    end
    for _,v in pairs(temp_w1) do
      if( not hash_w1[v]) then
        write1:insert(v)
        hash_w1[v] = true
      end
    end
    for _,v in pairs(temp_w2) do
      if( not hash_w2[v]) then
        write2:insert(v)
        hash_w2[v] = true
      end
    end
  end

  --Create the tasks
  local cell_pair_task, update_neighbours1 = generate_symmetric_pairwise_task_multikernel( kernels, read1, read2, write1, write2 )
  local cell_self_task, update_neighbours2 = generate_symmetric_self_task_multikernel( kernels, read1, read2, write1, write2 )
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
                --No periodic wrapping
                if (cell2_x >= 0 and cell2_x < x_count) 
                   and (cell2_y >= 0 and cell2_y < y_count)
                   and (cell2_z >= 0 and cell2_z < z_count) then
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
    [update_neighbours_quote];
end
return symmetric

end

-- ... should contain a list of symmetric kernel functions
function create_asymmetric_pairwise_runner_multikernel( config, cell_space, ... )
  local read1 = terralib.newlist()
  local read2 = terralib.newlist()
  local write1 = terralib.newlist()
  local write2 = terralib.newlist()
  if select("#", ...) == 0 then
    print("Attempting to make a 0 kernel task")
    return nil
  end
  local kernels = terralib.newlist()
  local hash_r1 = {}
  local hash_r2 = {}
  local hash_w1 = {}
  local hash_w2 = {}
  for i=1, select("#", ...) do
    kernels:insert(select(i, ...))
    local temp_r1, temp_r2, temp_w1, temp_w2 = compute_privileges.two_region_privileges(select(i, ...))
    --Merge the read/writes for this kernel with previous ones, keeping uniqueness
    for _,v in pairs(temp_r1) do
      if( not hash_r1[v]) then
        read1:insert(v)
        hash_r1[v] = true
      end
    end
    for _,v in pairs(temp_r2) do
      if( not hash_r2[v]) then
        read2:insert(v)
        hash_r2[v] = true
      end
    end
    for _,v in pairs(temp_w1) do
      if( not hash_w1[v]) then
        write1:insert(v)
        hash_w1[v] = true
      end
    end
    for _,v in pairs(temp_w2) do
      if( not hash_w2[v]) then
        write2:insert(v)
        hash_w2[v] = true
      end
    end
  end

  --Create the tasks
  local cell_pair_task, update_neighbours1 = generate_asymmetric_pairwise_task_multikernel( kernels, read1, read2, write1 )
  local cell_self_task, update_neighbours2 = generate_asymmetric_self_task_multikernel( kernels, read1, read2, write1)
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
                --No periodic wrapping
                if (cell2_x >= 0 and cell2_x < x_count) 
                   and (cell2_y >= 0 and cell2_y < y_count)
                   and (cell2_z >= 0 and cell2_z < z_count) then
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
    [update_neighbours_quote];
end
return asymmetric
end
-- ... should contain a list of per_particle kernel functions.
function run_per_particle_task_multikernel( config, cell_space, ...)
  local read1 = terralib.newlist()
  local read2 = terralib.newlist()
  local write1 = terralib.newlist()
  local write2 = terralib.newlist()
  if select("#", ...) == 0 then
    print("Attempting to make a 0 kernel task")
    return nil
  end
  local kernels = terralib.newlist()
  local hash_r1 = {}
  local hash_r2 = {}
  local hash_w1 = {}
  local hash_w2 = {}
  for i=1, select("#", ...) do
    kernels:insert(select(i, ...))
    local temp_r1, temp_r2, temp_w1, temp_w2 = compute_privileges.two_region_privileges(select(i, ...))
    --Merge the read/writes for this kernel with previous ones, keeping uniqueness 
    for _,v in pairs(temp_r1) do
      if( not hash_r1[v]) then
        read1:insert(v)
        hash_r1[v] = true
      end
    end
    for _,v in pairs(temp_r2) do
      if( not hash_r2[v]) then
        read2:insert(v)
        hash_r2[v] = true
      end
    end
    for _,v in pairs(temp_w1) do
      if( not hash_w1[v]) then
        write1:insert(v)
        hash_w1[v] = true
      end
    end
    for _,v in pairs(temp_w2) do
      if( not hash_w2[v]) then
        write2:insert(v)
        hash_w2[v] = true
      end
    end
  end
  local per_part_task, update_neighbours = generate_per_part_task_multikernel( kernels, read1, write1 )
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
--End of Multikernel Launchers-----------
-----------------------------------------

-----------------------------------------
--Invoke Framework-----------------------
-----------------------------------------

--By default, asking for PAIRWISE gives a SYMMETRIC_PAIRWISE operation
SYMMETRIC_PAIRWISE = 1
PAIRWISE = 1
ASYMMETRIC_PAIRWISE = 2
PER_PART = 3

BARRIER = 100
NO_BARRIER = 101

MULTI_KERNEL=1000
SINGLE_KERNEL=1001

local function is_safe_to_combine(kernel, combined_kernels)

  local pre_read1 = terralib.newlist()
  local pre_write1 = terralib.newlist()
  local hash_r1 = {}
  local hash_w1 = {}
  --Compute the read/write requirements for the already combined kernels
  for _, kernel in pairs(combined_kernels) do
    local temp_r1, temp_r2, temp_w1, temp_w2 = compute_privileges.two_region_privileges(kernel)
    --Merge the read/writes for this kernel with previous ones, keeping uniqueness
    for _,v in pairs(temp_r1) do
      if( not hash_r1[v]) then
        pre_read1:insert(v)
        hash_r1[v] = true
      end
    end
    for _,v in pairs(temp_r2) do
      if( not hash_r1[v]) then
        pre_read1:insert(v)
        hash_r1[v] = true
      end
    end
    for _,v in pairs(temp_w1) do
      if( not hash_w1[v]) then
        pre_write1:insert(v)
        hash_w1[v] = true
      end
    end
    for _,v in pairs(temp_w2) do
      if( not hash_w1[v]) then
        pre_write1:insert(v)
        hash_w1[v] = true
      end
    end
  end

local read1, read2, write1, write2 = compute_privileges.two_region_privileges( kernel )
local safe_to_combine = true
for _, v in pairs(write1) do
  if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" or v == "core_part_space.cutoff" then
    safe_to_combine = false
  end
  --Handle WaW or WaR dependencies
  if (hash_w1[v]) or (hash_r1[v]) then
    safe_to_combine = false
  end
end
for _, v in pairs(write2) do
  if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" or v == "core_part_space.cutoff" then
    safe_to_combine = false
  end
  --Handle WaW or WaR dependencies
  if (hash_w1[v]) or (hash_r1[v]) then
    safe_to_combine = false
  end
end
for _,v in pairs(read1) do
  --Handle RaW dependency
  if (hash_w1[v]) then
    safe_to_combine = false
  end
end
for _,v in pairs(read2) do
  --Handle RaW dependency
  if (hash_w1[v]) then
    safe_to_combine = false
  end
end
return safe_to_combine
end

function invoke_multikernel(config, ...)
  local end_barrier = true
  --Loop over the arguments and check for synchronicity (or other things to be added).
  for i= 1, select("#",...) do
    if select(i, ...) == NO_BARRIER then
      end_barrier = false
    end
  end
  local kernels = {}
  local last_type = -1
  local quote_list = terralib.newlist()

  --Loop through the inputs and find all the tables.
  for i=1, select("#", ...) do
    local v = select(i, ...)
    if type(v) == "table" then
      if v[1] == nil or v[2] == nil then
        print("The arguments to invoke were incorrect")
        os.exit(1)
      end
      local func = v[1]
      local type_iterate = v[2]
 --I think we can refactor this using some functions to make the code cleaner, and just check if last_type == type_iterate. AC
      if type_iterate == SYMMETRIC_PAIRWISE then
        local safe_to_combine, can_be_combined = is_safe_to_combine(func, kernels)
        if safe_to_combine and last_type == SYMMETRIC_PAIRWISE then
          table.insert(kernels, func)
        elseif safe_to_combine then
          if #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
            quote_list:insert( create_asymmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == PER_PART then
            quote_list:insert( run_per_particle_task_multikernel( config, neighbour_init.cell_partition, unpack(kernels) ) )
          end

          kernels = {}
          table.insert(kernels, func)
          last_type = SYMMETRIC_PAIRWISE
        else --GENERATE PREVIOUS AND CURRENT
          if #kernels > 0 and last_type == SYMMETRIC_PAIRWISE then
            quote_list:insert( create_symmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
            quote_list:insert( create_asymmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == PER_PART then
            quote_list:insert( run_per_particle_task_multikernel( config, neighbour_init.cell_partition, unpack(kernels) ) )
          end
          if can_be_combined then
            kernels = {}
            table.insert(kernels, func)
            last_type = SYMMETRIC_PAIRWISE
          else
            quote_list:insert( create_symmetric_pairwise_runner(func, config, neighbour_init.cell_partition))
            kernels = {}
            last_type = -1
          end
        end
      elseif type_iterate == ASYMMETRIC_PAIRWISE then
        local safe_to_combine, can_be_combined = is_safe_to_combine(func, kernels)
        if safe_to_combine and last_type == ASYMMETRIC_PAIRWISE then
          table.insert(kernels, func)
        elseif safe_to_combine then
          if #kernels > 0 and last_type == SYMMETRIC_PAIRWISE then
            quote_list:insert( create_symmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == PER_PART then
            quote_list:insert( run_per_particle_task_multikernel( config, neighbour_init.cell_partition, unpack(kernels) ) )
          end

          kernels = {}
          table.insert(kernels, func)
          last_type = ASYMMETRIC_PAIRWISE
        else --GENERATE PREVIOUS AND CURRENT
          if #kernels > 0 and last_type == SYMMETRIC_PAIRWISE then
            quote_list:insert( create_symmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
            quote_list:insert( create_asymmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == PER_PART then
            quote_list:insert( run_per_particle_task_multikernel( config, neighbour_init.cell_partition, unpack(kernels) ) )
          end

          if can_be_combined then
            kernels = {}
            table.insert(kernels, func)
            last_type = ASYMMETRIC_PAIRWISE
          else
            quote_list:insert( create_asymmetric_pairwise_runner(func, config, neighbour_init.cell_partition))
            kernels = {}
            last_type = -1
          end
        end
      elseif type_iterate == PER_PART then
        local safe_to_combine, can_be_combined = is_safe_to_combine(func, kernels)
        if safe_to_combine and last_type == PER_PART then
          table.insert(kernels, func)
        elseif safe_to_combine then
          if #kernels > 0 and last_type == SYMMETRIC_PAIRWISE then
            quote_list:insert( create_symmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
            quote_list:insert( create_asymmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          end

          kernels = {}
          table.insert(kernels, func)
          last_type = PER_PART
        else --GENERATE PREVIOUS AND CURRENT
          if #kernels > 0 and last_type == SYMMETRIC_PAIRWISE then
            quote_list:insert( create_symmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
            quote_list:insert( create_asymmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == PER_PART then
            quote_list:insert( run_per_particle_task_multikernel( config, neighbour_init.cell_partition, unpack(kernels) ) )
          end
          if can_be_combined then
            kernels = {}
            table.insert(kernels, func)
            last_type = PER_PART
          else
            quote_list:insert(  run_per_particle_task( func, config, neighbour_init.cell_partition ) )
            kernels = {}
            last_type = -1
          end
        end
      else
        print("The kernel type passed to invoke was not recognized")
        os.exit(1)
      end
    end
  end
  --Generate the final quote
  if #kernels > 0 and last_type == SYMMETRIC_PAIRWISE then
    quote_list:insert( create_symmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
  elseif #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
    quote_list:insert( create_asymmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
  elseif #kernels > 0 and last_type == PER_PART then
    quote_list:insert( run_per_particle_task_multikernel( config, neighbour_init.cell_partition, unpack(kernels) ) )
  end
  local barrier_quote = rquote
  end
  if end_barrier then
    barrier_quote = rquote
      c.legion_runtime_issue_execution_fence(__runtime(), __context())
    end
  end
  local invoke_quote = rquote
    [quote_list];
    [barrier_quote];
  end
  return invoke_quote
end

function invoke_per_kernel(config, ...)
  local end_barrier = true
  --Loop over the arguments and check for synchronicity (or other things to be added).
  for i= 1, select("#",...) do
    if select(i, ...) == NO_BARRIER then
      end_barrier = false
    end
  end
  
  local quote_list = terralib.newlist()
  --Loop through the inputs and find all the tables. 
  for i= 1, select("#",...) do
     local v = select(i, ...)
    if type(v) == "table" then
      if v[1] == nil or v[2] == nil then
        print("The arguments to invoke were incorrect")
        os.exit(1)
      end
      local func = v[1]
      local type_iterate = v[2]
      if type_iterate == SYMMETRIC_PAIRWISE then
        quote_list:insert( create_symmetric_pairwise_runner(func, config, neighbour_init.cell_partition) )
      elseif type_iterate == ASYMMETRIC_PAIRWISE then
        quote_list:insert( create_asymmetric_pairwise_runner(func, config, neighbour_init.cell_partition) )
      elseif type_iterate == PER_PART then
        quote_list:insert( run_per_particle_task(func, config,  neighbour_init.cell_partition) ) 
      else
        print("The kernel type passed to invoke was not recognized")
        os.exit(1)
      end
    end
  end
  local barrier_quote = rquote

  end
  if end_barrier then
    barrier_quote = rquote
      c.legion_runtime_issue_execution_fence(__runtime(), __context())
    end
  end 
  local invoke_quote = rquote
    [quote_list];
    [barrier_quote];
  end
return invoke_quote
end


function invoke(config, ...)
  local multi_kernel = true
  --Loop over the arguments and check for synchronicity (or other things to be added).
  for i= 1, select("#",...) do
    if select(i, ...) == SINGLE_KERNEL then
      multi_kernel = false
    end
  end
  if multi_kernel then
    return invoke_multikernel(config, ...)
  else
    return invoke_per_kernel(config, ...)
  end
end

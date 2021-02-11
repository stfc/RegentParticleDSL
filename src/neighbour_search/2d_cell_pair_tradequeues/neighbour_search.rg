-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

require("src/neighbour_search/2d_cell_pair_tradequeues/cell")
local compute_privileges = require("src/utils/compute_privilege")
local format = require("std/format")
local string_to_field_path = require("src/utils/string_to_fieldpath")
local coherence_compute = require("src/utils/coherence_computing")
local privilege_lists = require("src/utils/privilege_lists")
local c = regentlib.c

local abs = regentlib.fabs(double)
local ceil = regentlib.ceil(double)


--This function assumes the cutoff is the same for both particles
function generate_symmetric_pairwise_task_multikernel( kernel_list, read1, read2, write1, write2, reduc1, reduc2 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local update_neighbours, read1_privs, read2_privs, write1_privs, write2_privs, reduc1_privs, reduc2_privs =
      privilege_lists.get_privileges_symmetric_pair_task( parts1, parts2, read1, read2,
              write1, write2, reduc1, reduc2 )
local coherences = coherence_compute.compute_coherences_pair_task(update_neighbours, parts1, parts2)

local __demand(__leaf) task pairwise_task([parts1], [parts2], config : region(ispace(int1d), config_type))
  where [read1_privs], [read2_privs], [write1_privs], [write2_privs], [reduc1_privs], [reduc2_privs],
  reads(config), reads(parts1.core_part_space.{pos_x, pos_y, cutoff}),
  reads(parts2.core_part_space.{pos_x, pos_y, cutoff}), reads(parts1.neighbour_part_space._valid), reads(parts2.neighbour_part_space._valid),
 [coherences] do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   --For the tradequeue implementation, we only care about the validity of the particles
   for part1 in [parts1].ispace do
     if [parts1][part1].neighbour_part_space._valid then
       for part2 in [parts2].ispace do
         if [parts2][part2].neighbour_part_space._valid then
           --Compute particle distance
           var dx = [parts1][part1].core_part_space.pos_x - [parts2][part2].core_part_space.pos_x
           var dy = [parts1][part1].core_part_space.pos_y - [parts2][part2].core_part_space.pos_y
           if (dx > half_box_x) then dx = dx - box_x end
           if (dy > half_box_y) then dy = dy - box_y end
           if (dx <-half_box_x) then dx = dx + box_x end
           if (dy <-half_box_y) then dy = dy + box_y end
           var cutoff2 = regentlib.fmax([parts1][part1].core_part_space.cutoff, [parts2][part2].core_part_space.cutoff)
           cutoff2 = cutoff2 * cutoff2
           var r2 = dx*dx + dy*dy
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
function generate_symmetric_pairwise_task( kernel_name, read1, read2, write1, write2, reduc1, reduc2 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local update_neighbours, read1_privs, read2_privs, write1_privs, write2_privs, reduc1_privs, reduc2_privs =
      privilege_lists.get_privileges_symmetric_pair_task( parts1, parts2, read1, read2,
              write1, write2, reduc1, reduc2 )
local coherences = coherence_compute.compute_coherences_pair_task(update_neighbours, parts1, parts2)

local __demand(__leaf) task pairwise_task([parts1], [parts2], config : region(ispace(int1d), config_type))
  where [read1_privs], [read2_privs], [write1_privs], [write2_privs], [reduc1_privs], [reduc2_privs],
  reads(config), reads(parts1.core_part_space.{pos_x, pos_y, cutoff}),
  reads(parts2.core_part_space.{pos_x, pos_y, cutoff}), reads(parts1.neighbour_part_space._valid), reads(parts2.neighbour_part_space._valid),
 [coherences] do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   --For the tradequeue implementation, we only care about the validity of the particles
   for part1 in [parts1].ispace do
     if [parts1][part1].neighbour_part_space._valid then
       for part2 in [parts2].ispace do
         if [parts2][part2].neighbour_part_space._valid then
           --Compute particle distance
           var dx = [parts1][part1].core_part_space.pos_x - [parts2][part2].core_part_space.pos_x
           var dy = [parts1][part1].core_part_space.pos_y - [parts2][part2].core_part_space.pos_y
           if (dx > half_box_x) then dx = dx - box_x end
           if (dy > half_box_y) then dy = dy - box_y end
           if (dx <-half_box_x) then dx = dx + box_x end
           if (dy <-half_box_y) then dy = dy + box_y end
           var cutoff2 = regentlib.fmax([parts1][part1].core_part_space.cutoff, [parts2][part2].core_part_space.cutoff)
           cutoff2 = cutoff2 * cutoff2
           var r2 = dx*dx + dy*dy
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

function generate_asymmetric_pairwise_task_multikernel( kernel_list, read1, read2, write1, reduc1 )
--Asymmetric kernel can only write to part1
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local update_neighbours, read1_privs, read2_privs, write1_privs, reduc1_privs =
      privilege_lists.get_privileges_asymmetric_pair_task( parts1, parts2, read1, read2,
              write1, reduc1 )
local coherences = coherence_compute.compute_coherences_pair_task(update_neighbours, parts1, parts2)

local __demand(__leaf) task pairwise_task([parts1], [parts2],  config : region(ispace(int1d), config_type))
  where [read1_privs], [read2_privs], [write1_privs], [reduc1_privs],
  reads(config), reads(parts1.core_part_space.{pos_x, pos_y, cutoff}),
  reads(parts2.core_part_space.{pos_x, pos_y, cutoff}), reads(parts1.neighbour_part_space._valid), reads(parts2.neighbour_part_space._valid),
  [coherences] do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   for part1 in [parts1].ispace do
     if [parts1][part1].neighbour_part_space._valid then
       for part2 in [parts2].ispace do
         if [parts2][part2].neighbour_part_space._valid then
           --Compute particle distance
           var dx = [parts1][part1].core_part_space.pos_x - [parts2][part2].core_part_space.pos_x
           var dy = [parts1][part1].core_part_space.pos_y - [parts2][part2].core_part_space.pos_y
           if (dx > half_box_x) then dx = dx - box_x end
           if (dy > half_box_y) then dy = dy - box_y end
           if (dx <-half_box_x) then dx = dx + box_x end
           if (dy <-half_box_y) then dy = dy + box_y end
           var cutoff2 = [parts1][part1].core_part_space.cutoff
           cutoff2 = cutoff2 * cutoff2
           var r2 = dx*dx + dy*dy
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

function generate_asymmetric_pairwise_task( kernel_name, read1, read2, write1, reduc1 )
--Asymmetric kernel can only write to part1
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local update_neighbours, read1_privs, read2_privs, write1_privs, reduc1_privs =
      privilege_lists.get_privileges_asymmetric_pair_task( parts1, parts2, read1, read2,
              write1, reduc1 )
local coherences = coherence_compute.compute_coherences_pair_task(update_neighbours, parts1, parts2)

local __demand(__leaf) task pairwise_task([parts1], [parts2],  config : region(ispace(int1d), config_type))
  where [read1_privs], [read2_privs], [write1_privs], [reduc1_privs], reads(config), reads(parts1.core_part_space.{pos_x, pos_y, cutoff}),
  reads(parts2.core_part_space.{pos_x, pos_y, cutoff}), reads(parts1.neighbour_part_space._valid), reads(parts2.neighbour_part_space._valid),
  [coherences] do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   for part1 in [parts1].ispace do
     if [parts1][part1].neighbour_part_space._valid then
       for part2 in [parts2].ispace do
         if [parts2][part2].neighbour_part_space._valid then
           --Compute particle distance
           var dx = [parts1][part1].core_part_space.pos_x - [parts2][part2].core_part_space.pos_x
           var dy = [parts1][part1].core_part_space.pos_y - [parts2][part2].core_part_space.pos_y
           if (dx > half_box_x) then dx = dx - box_x end
           if (dy > half_box_y) then dy = dy - box_y end
           if (dx <-half_box_x) then dx = dx + box_x end
           if (dy <-half_box_y) then dy = dy + box_y end
           var cutoff2 = [parts1][part1].core_part_space.cutoff
           cutoff2 = cutoff2 * cutoff2
           var r2 = dx*dx + dy*dy
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

function generate_symmetric_self_task_multikernel( kernel_list, read1, read2, write1, write2, reduc1, reduc2 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local update_neighbours, read1_privs, write1_privs, reduc1_privs = privilege_lists.get_privileges_self_task( parts1, read1, read2, write1, write2, reduc1, reduc2 )
local coherences = coherence_compute.compute_coherences_self_task(update_neighbours, parts1)

local __demand(__leaf) task self_task([parts1], config : region(ispace(int1d),config_type)) where
  [read1_privs], [write1_privs], [reduc1_privs], reads(parts1.core_part_space.{pos_x, pos_y, cutoff}),
                                 reads(parts1.neighbour_part_space._valid), reads(config),
   [coherences] do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   for part1 in [parts1].ispace do
     if [parts1][part1].neighbour_part_space._valid then
       for part2 in [parts1].ispace do
         if [parts1][part2].neighbour_part_space._valid then
           --Compute particle distance
           if(part1 < part2) then
             var dx = [parts1][part1].core_part_space.pos_x - [parts1][part2].core_part_space.pos_x
             var dy = [parts1][part1].core_part_space.pos_y - [parts1][part2].core_part_space.pos_y
             if (dx > half_box_x) then dx = dx - box_x end
             if (dy > half_box_y) then dy = dy - box_y end
             if (dx <-half_box_x) then dx = dx + box_x end
             if (dy <-half_box_y) then dy = dy + box_y end
             var cutoff2 = regentlib.fmax([parts1][part1].core_part_space.cutoff, [parts1][part2].core_part_space.cutoff)
             cutoff2 = cutoff2 * cutoff2
             var r2 = dx*dx + dy*dy
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
function generate_symmetric_self_task( kernel_name, read1, read2, write1, write2, reduc1, reduc2 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local update_neighbours, read1_privs, write1_privs, reduc1_privs = privilege_lists.get_privileges_self_task( parts1, read1, read2, write1, write2, reduc1, reduc2 )
local coherences = coherence_compute.compute_coherences_self_task(update_neighbours, parts1)

local __demand(__leaf) task self_task([parts1], config : region(ispace(int1d),config_type)) where
  [read1_privs], [write1_privs], [reduc1_privs], reads(parts1.core_part_space.{pos_x, pos_y, cutoff}),
                                 reads(parts1.neighbour_part_space._valid), reads(config),
   [coherences] do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   for part1 in [parts1].ispace do
     if [parts1][part1].neighbour_part_space._valid then
       for part2 in [parts1].ispace do
         if [parts1][part2].neighbour_part_space._valid then
           --Compute particle distance
           if(part1 < part2) then
             var dx = [parts1][part1].core_part_space.pos_x - [parts1][part2].core_part_space.pos_x
             var dy = [parts1][part1].core_part_space.pos_y - [parts1][part2].core_part_space.pos_y
             if (dx > half_box_x) then dx = dx - box_x end
             if (dy > half_box_y) then dy = dy - box_y end
             if (dx <-half_box_x) then dx = dx + box_x end
             if (dy <-half_box_y) then dy = dy + box_y end
             var cutoff2 = regentlib.fmax([parts1][part1].core_part_space.cutoff, [parts1][part2].core_part_space.cutoff)
             cutoff2 = cutoff2 * cutoff2
             var r2 = dx*dx + dy*dy
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

function generate_asymmetric_self_task_multikernel( kernel_list, read1, read2, write1, reduc1 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local update_neighbours, read1_privs, write1_privs, reduc1_privs = privilege_lists.get_privileges_self_task( parts1, read1, read2, write1, {}, reduc1, {} )
local coherences = coherence_compute.compute_coherences_self_task(update_neighbours, parts1)

local __demand(__leaf) task self_task([parts1], config : region(ispace(int1d), config_type)) where
   [read1_privs], [write1_privs], [reduc1_privs], reads(parts1.core_part_space.{pos_x, pos_y, cutoff}),
                                  reads(parts1.neighbour_part_space._valid), reads(config),
   [coherences] do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   for part1 in [parts1].ispace do
     if [parts1][part1].neighbour_part_space._valid then
       for part2 in [parts1].ispace do
         if [parts1][part2].neighbour_part_space._valid then
           --Compute particle distance
           if(part1 ~= part2) then
             var dx = [parts1][part1].core_part_space.pos_x - [parts1][part2].core_part_space.pos_x
             var dy = [parts1][part1].core_part_space.pos_y - [parts1][part2].core_part_space.pos_y
             if (dx > half_box_x) then dx = dx - box_x end
             if (dy > half_box_y) then dy = dy - box_y end
             if (dx <-half_box_x) then dx = dx + box_x end
             if (dy <-half_box_y) then dy = dy + box_y end
             var cutoff2 = [parts1][part1].core_part_space.cutoff
             cutoff2 = cutoff2 * cutoff2
             var r2 = dx*dx + dy*dy
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
function generate_asymmetric_self_task( kernel_name, read1, read2, write1, reduc1 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local update_neighbours, read1_privs, write1_privs, reduc1_privs = privilege_lists.get_privileges_self_task( parts1, read1, read2, write1, {}, reduc1, {} )
local coherences = coherence_compute.compute_coherences_self_task(update_neighbours, parts1)

local __demand(__leaf) task self_task([parts1], config : region(ispace(int1d), config_type)) where
   [read1_privs], [write1_privs], [reduc1_privs], reads(parts1.core_part_space.{pos_x, pos_y, cutoff}), 
                                  reads(parts1.neighbour_part_space._valid), reads(config),
   [coherences] do
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   for part1 in [parts1].ispace do
     if [parts1][part1].neighbour_part_space._valid then
       for part2 in [parts1].ispace do
         if [parts1][part2].neighbour_part_space._valid then
           --Compute particle distance
           if(part1 ~= part2) then
             var dx = [parts1][part1].core_part_space.pos_x - [parts1][part2].core_part_space.pos_x
             var dy = [parts1][part1].core_part_space.pos_y - [parts1][part2].core_part_space.pos_y
             if (dx > half_box_x) then dx = dx - box_x end
             if (dy > half_box_y) then dy = dy - box_y end
             if (dx <-half_box_x) then dx = dx + box_x end
             if (dy <-half_box_y) then dy = dy + box_y end
             var cutoff2 = [parts1][part1].core_part_space.cutoff
             cutoff2 = cutoff2 * cutoff2
             var r2 = dx*dx + dy*dy
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
task cell_greater_equal(cell1 : int2d, cell2 : int2d) : bool
  var returnval : bool = false
  if (cell1.x > cell2.x or (cell1.x == cell2.x and cell1.y > cell2.y) ) then
    returnval =  true
  end
  return returnval
end

function create_symmetric_pairwise_runner( kernel_name, config, cell_space )
local read1, read2, write1, write2, reduc1, reduc2 = compute_privileges.two_region_privileges( kernel_name )
local cell_pair_task, update_neighbours1 = generate_symmetric_pairwise_task( kernel_name, read1, read2, write1, write2, reduc1, reduc2 )
local cell_self_task, update_neighbours2 = generate_symmetric_self_task( kernel_name, read1, read2, write1, write2, reduc1, reduc2 )
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

    --Compute cell radii
    var cutoff = config[0].neighbour_config.max_cutoff
    var x_radii : int = ceil( cutoff / config[0].neighbour_config.cell_dim_x )
    var y_radii : int = ceil( cutoff / config[0].neighbour_config.cell_dim_y )
    --We want to trace all the tasks here
    __demand(__trace)
    do
        --Spawn all the self tasks first, which can be an index launch!
        __demand(__index_launch)
        for cell1 in cell_space.colors do
            cell_self_task(cell_space[cell1], config)
        end
        --Now spawn all the pairs
        for cell1 in cell_space.colors do
            --Loops non inclusive, positive only direction.
            for x = -x_radii, x_radii+1 do
              for y = -y_radii, y_radii+1 do
                if(not (x == 0 and y == 0) ) then
                  var cell2_x = cell1.x + x
                  var cell2_y = cell1.y + y
                  if cell2_x < 0 then
                    cell2_x = cell2_x + x_count
                  end
                  if cell2_y < 0 then
                    cell2_y = cell2_y + y_count
                  end
                  var cell2 : int2d = int2d({ (cell2_x)%x_count, (cell2_y)%y_count })
                  --Weird if statement to handle max_cutoff >= half the boxsize
                  if( cell_greater_equal(cell1, cell2)
                  or ((x < 0 and cell2.x <= cell1.x + x_radii and cell2.x > cell1.x) or
                                                          (y < 0 and cell2.y <= cell1.y + y_radii and cell2.y > cell1.y))
                  or ((x > 0 and cell2.x >= cell1.x - x_radii and cell2.x < cell1.x) or (y > 0 and cell2.y >= cell1.y - y_radii and cell2.y < cell1.y))
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
local read1, read2, write1, write2, reduc1, reduc2 = compute_privileges.two_region_privileges( kernel_name )
--While write2 is computed, asymmetric kernels are not allowed to write to write2
local cell_pair_task, update_neighbours1 = generate_asymmetric_pairwise_task( kernel_name, read1, read2, write1, reduc1, reduc2 )
local cell_self_task, update_neighbours2 = generate_asymmetric_self_task( kernel_name, read1, read2, write1, reduc1, reduc2 )
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

    --Compute cell radii
    var cutoff = config[0].neighbour_config.max_cutoff
    var x_radii : int = ceil( cutoff / config[0].neighbour_config.cell_dim_x )
    var y_radii : int = ceil( cutoff / config[0].neighbour_config.cell_dim_y )
    __demand(__trace)
    do
        --Launch self tasks then pair tasks
        __demand(__index_launch)
        for cell1 in cell_space.colors do
            cell_self_task(cell_space[cell1], config)
        end
        for cell1 in cell_space.colors do
            --Loops non inclusive, positive only direction.
            for x = -x_radii, x_radii+1 do
              for y = -y_radii, y_radii+1 do
                if(not (x == 0 and y == 0 ) ) then
                  var cell2_x = cell1.x + x
                  var cell2_y = cell1.y + y
                  if cell2_x < 0 then
                    cell2_x = cell2_x + x_count
                  end
                  if cell2_y < 0 then
                    cell2_y = cell2_y + y_count
                  end
                  var cell2 : int2d = int2d({ (cell2_x)%x_count, (cell2_y)%y_count })
                  --if statement to handle max_cutoff >= half the boxsize
                      --Handle radius overlap
                  if( cell_greater_equal(cell1, cell2)
                  or ((x < 0 and cell2.x <= cell1.x + x_radii and cell2.x > cell1.x) or
                                                          (y < 0 and cell2.y <= cell1.y + y_radii and cell2.y > cell1.y))
                  or ((x > 0 and cell2.x >= cell1.x - x_radii and cell2.x < cell1.x) or (y > 0 and cell2.y >= cell1.y - y_radii and cell2.y < cell1.y))
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

--Generate a task to be executed on every particle in the system, using a list of kernels
function generate_per_part_task_multikernel( kernel_list, read1, write1, reduc1 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local update_neighbours, read1_privs, write1_privs, reduc1_privs = privilege_lists.get_privileges_self_task( parts1, read1, {}, write1, {}, reduc1, {} )
local coherences = coherence_compute.compute_coherences_self_task(update_neighbours, parts1)

local __demand(__leaf) task pairwise_task([parts1], config : region(ispace(int1d), config_type)) where
   [read1_privs], [write1_privs], [reduc1_privs], reads(config), [coherences],
   reads( parts1.neighbour_part_space._valid ) 
   do
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
function generate_per_part_task( kernel_name, read1, write1, reduc1 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local update_neighbours, read1_privs, write1_privs, reduc1_privs = privilege_lists.get_privileges_self_task( parts1, read1, {}, write1, {}, reduc1, {} )
local coherences = coherence_compute.compute_coherences_self_task(update_neighbours, parts1)

local __demand(__leaf) task per_part_task([parts1], config : region(ispace(int1d), config_type)) where
   [read1_privs], [write1_privs], [reduc1_privs], reads(config), [coherences],
   reads( parts1.neighbour_part_space._valid )
   do
   for part1 in [parts1].ispace do
     if [parts1][part1].neighbour_part_space._valid then
         [kernel_name(rexpr [parts1][part1] end, rexpr config[0] end)]
     end
   end
end
return per_part_task, update_neighbours
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

local read1, read2, write1, write2, reduc1, reduc2 = compute_privileges.two_region_privileges(kernel_name)
local need_tradequeues = false
local per_part_task, update_neighbours = generate_per_part_task( kernel_name, read1, write1, reduc1 )
local update_neighbours_quote = rquote

end
if update_neighbours then
  local temp_variables = {}
  temp_variables.config = config
  update_neighbours_quote = neighbour_init.update_cells(temp_variables)
end
local runner = rquote

    --For each cell, call the task!
    __demand(__index_launch)
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
  local reduc1 = terralib.newlist()
  local reduc2 = terralib.newlist()
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
    local temp_r1, temp_r2, temp_w1, temp_w2, temp_re1, temp_re2 = compute_privileges.two_region_privileges(select(i, ...))
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
    for _, v in pairs(temp_re1) do
      local exists = false
      for _, v2 in pairs(reduc1) do
        if v[1] == v2[1] and v[2] == v2[2] then
          exists = true
        end
      end
      if not exists then
        reduc1:insert(v)
      end
    end
    for _, v in pairs(temp_re2) do
      local exists = false
      for _, v2 in pairs(reduc2) do
        if v[1] == v2[1] and v[2] == v2[2] then
          exists = true
        end
      end
      if not exists then
        reduc2:insert(v)
      end
    end
  end

  --Create the tasks
  local cell_pair_task, update_neighbours1 = generate_symmetric_pairwise_task_multikernel( kernels, read1, read2, write1, write2, reduc1, reduc2 )
  local cell_self_task, update_neighbours2 = generate_symmetric_self_task_multikernel( kernels, read1, read2, write1, write2, reduc1, reduc2 )
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

    --Compute cell radii
    var cutoff = config[0].neighbour_config.max_cutoff
    var x_radii : int = ceil( cutoff / config[0].neighbour_config.cell_dim_x )
    var y_radii : int = ceil( cutoff / config[0].neighbour_config.cell_dim_y )
    for cell1 in cell_space.colors do
        cell_self_task(cell_space[cell1], config)
        --Loops non inclusive, positive only direction.
        for x = -x_radii, x_radii+1 do
          for y = -y_radii, y_radii+1 do
            if(not (x == 0 and y == 0) ) then
              var cell2_x = cell1.x + x
              var cell2_y = cell1.y + y
              if cell2_x < 0 then
                cell2_x = cell2_x + x_count
              end
              if cell2_y < 0 then
                cell2_y = cell2_y + y_count
              end
              var cell2 : int2d = int2d({ (cell2_x)%x_count, (cell2_y)%y_count })
              --Weird if statement to handle max_cutoff >= half the boxsize
              if( cell_greater_equal(cell1, cell2)
              or ((x < 0 and cell2.x <= cell1.x + x_radii and cell2.x > cell1.x) or
                                                      (y < 0 and cell2.y <= cell1.y + y_radii and cell2.y > cell1.y))
              or ((x > 0 and cell2.x >= cell1.x - x_radii and cell2.x < cell1.x) or (y > 0 and cell2.y >= cell1.y - y_radii and cell2.y < cell1.y))
              ) then
              else
                --symmetric
                cell_pair_task(cell_space[cell1], cell_space[cell2], config)
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
  local reduc1 = terralib.newlist()
  local reduc2 = terralib.newlist()
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
    local temp_r1, temp_r2, temp_w1, temp_w2, temp_re1, temp_re2 = compute_privileges.two_region_privileges(select(i, ...))
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
    for _, v in pairs(temp_re1) do
      local exists = false                                                                                                                                                       for _, v2 in pairs(reduc1) do
        if v[1] == v2[1] and v[2] == v2[2] then
          exists = true
        end
      end
      if not exists then
        reduc1:insert(v)
      end
    end
  end

  --Create the tasks
  local cell_pair_task, update_neighbours1 = generate_asymmetric_pairwise_task_multikernel( kernels, read1, read2, write1, temp_re1 )
  local cell_self_task, update_neighbours2 = generate_asymmetric_self_task_multikernel( kernels, read1, read2, write1, temp_re1 )
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

    --Compute cell radii
    var cutoff = config[0].neighbour_config.max_cutoff
    var x_radii : int = ceil( cutoff / config[0].neighbour_config.cell_dim_x )
    var y_radii : int = ceil( cutoff / config[0].neighbour_config.cell_dim_y )
    for cell1 in cell_space.colors do
        cell_self_task(cell_space[cell1], config)
        --Loops non inclusive, positive only direction.
        for x = -x_radii, x_radii+1 do
          for y = -y_radii, y_radii+1 do
            if(not (x == 0 and y == 0) ) then
              var cell2_x = cell1.x + x
              var cell2_y = cell1.y + y
              if cell2_x < 0 then
                cell2_x = cell2_x + x_count
              end
              if cell2_y < 0 then
                cell2_y = cell2_y + y_count
              end
              var cell2 : int2d = int2d({ (cell2_x)%x_count, (cell2_y)%y_count })
              --if statement to handle max_cutoff >= half the boxsize
                  --Handle radius overlap
              if( cell_greater_equal(cell1, cell2)
              or ((x < 0 and cell2.x <= cell1.x + x_radii and cell2.x > cell1.x) or
                                                      (y < 0 and cell2.y <= cell1.y + y_radii and cell2.y > cell1.y)) 
              or ((x > 0 and cell2.x >= cell1.x - x_radii and cell2.x < cell1.x) or (y > 0 and cell2.y >= cell1.y - y_radii and cell2.y < cell1.y))
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
  local reduc1 = terralib.newlist()
  local reduc2 = terralib.newlist()
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
    local temp_r1, temp_r2, temp_w1, temp_w2, temp_re1, temp_re2 = compute_privileges.two_region_privileges(select(i, ...))
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
    for _, v in pairs(temp_re1) do
      local exists = false
      for _, v2 in pairs(reduc1) do
        if v[1] == v2[1] and v[2] == v2[2] then
          exists = true
        end
      end
      if not exists then
        reduc1:insert(v)
      end
    end
    for _, v in pairs(temp_re2) do
      local exists = false
      for _, v2 in pairs(reduc2) do
        if v[1] == v2[1] and v[2] == v2[2] then
          exists = true
        end
      end
      if not exists then
        reduc2:insert(v)
      end
    end
  end
  local per_part_task, update_neighbours = generate_per_part_task_multikernel( kernels, read1, write1, reduc1 )
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

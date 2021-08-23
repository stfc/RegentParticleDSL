-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

require("src/neighbour_search/HP_cell_pair_tradequeues/cell")
local compute_privileges = require("src/utils/compute_privilege")
local format = require("std/format")
local string_to_field_path = require("src/utils/string_to_fieldpath")
local coherence_compute = require("src/utils/coherence_computing")
local privilege_lists = require("src/utils/privilege_lists")
local c = regentlib.c

local abs = regentlib.fabs(double)
local ceil = regentlib.ceil(double)
local sqrt = regentlib.sqrt(double)
---------------------------------------------------------------------------------------------
----------------------------- Utility function to wrap indices ------------------------------
---------------------------------------------------------------------------------------------
local terra wrap( val : int, size : int) : int
    if val < 0 then
        val = val + size
    end
    val = val % size
    return val
end
----------------------------------------------------------------------------------------------



--This function assumes the cutoff is the same for both particles
function generate_symmetric_pairwise_task( kernel_name, read1, read2, read3, write1, write2, reduc1, reduc2, reduc3 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local config = regentlib.newsymbol(region(ispace(int1d), config_type), "config")
local update_neighbours, read1_privs, read2_privs, read3_privs, write1_privs, write2_privs, reduc1_privs, reduc2_privs, reduc3_privs = 
      privilege_lists.get_privileges_symmetric_pair_task( parts1, parts2, config, read1, read2, read3,
              write1, write2, reduc1, reduc2, reduc3 )
--Also need to include the space in the reads part of the config
--read3_privs:insert( regentlib.privilege(regentlib.reads, config, string_to_field_path.get_field_path("space")))
local coherences = coherence_compute.compute_coherences_pair_task(update_neighbours, parts1, parts2)
local __demand(__leaf) task pairwise_task([parts1], [parts2], [config] )
  where [read1_privs], [read2_privs], [read3_privs], [write1_privs], [write2_privs], [reduc1_privs], [reduc2_privs], [reduc3_privs], reads(config.space), 
  reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}),
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
             [kernel_name(rexpr [parts1][part1] end, rexpr [parts2][part2] end, rexpr r2 end, rexpr config[0] end)];
--             [kernel_list:map( function(kernel)
--               return kernel(rexpr [parts1][part1] end, rexpr [parts2][part2] end, rexpr r2 end)
--             end)];
           end
         end
       end
     end
   end
end
return pairwise_task, update_neighbours
end

function generate_asymmetric_pairwise_task( kernel_name, read1, read2, read3, write1, reduc1, reduc3 )
--Asymmetric kernel can only write to part1
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local config = regentlib.newsymbol(region(ispace(int1d), config_type), "config")
local  update_neighbours, read1_privs, read2_privs, read3_privs, write1_privs, reduc1_privs, reduc3_privs =
      privilege_lists.get_privileges_asymmetric_pair_task( parts1, parts2, config, read1, read2, read3,
              write1, reduc1, reduc3 )
local coherences = coherence_compute.compute_coherences_pair_task(update_neighbours, parts1, parts2)
local __demand(__leaf) task asym_pairwise_task([parts1], [parts2], [config] )
  where [read1_privs], [read2_privs], [read3_privs], [write1_privs], [reduc1_privs], [reduc3_privs], reads(config.space), reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}),
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
             [kernel_name(rexpr [parts1][part1] end, rexpr [parts2][part2] end, rexpr r2 end, rexpr config[0] end)]
           end
         end
       end
     end
   end
end
return asym_pairwise_task, update_neighbours
end

-----------------------------------------
--End of Pair Tasks ---------------------
-----------------------------------------


--Generate a self task
--This function assumes the cutoff is the same for both particles
function generate_symmetric_self_task( kernel_name, read1, read2, read3, write1, write2, reduc1, reduc2, reduc3 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local config = regentlib.newsymbol(region(ispace(int1d), config_type), "config")
local update_neighbours, read1_privs, readconf_privs, write1_privs, reduc1_privs, reducconf_privs = privilege_lists.get_privileges_self_task(  parts1, config, read1, read2, read3, write1, write2, reduc1, reduc2, reduc3)
local coherences = coherence_compute.compute_coherences_self_task(update_neighbours, parts1)

local __demand(__leaf) task self_task([parts1], [config]) where
  [read1_privs], [readconf_privs], [write1_privs], [reduc1_privs], [reducconf_privs], reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}),
                                 reads(parts1.neighbour_part_space._valid), reads(config.space),
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
               [kernel_name(rexpr [parts1][part1] end, rexpr [parts1][part2] end, rexpr r2 end, rexpr config[0] end)]
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
function generate_asymmetric_self_task( kernel_name, read1, read2, read3, write1, reduc1, reduc3 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local config = regentlib.newsymbol(region(ispace(int1d), config_type), "config")
local update_neighbours, read1_privs, readconf_privs, write1_privs, reduc1_privs, reducconf_privs = 
                                                                   privilege_lists.get_privileges_self_task( parts1, config, read1, read2, read3,
                                                                                                             write1, {}, reduc1, {}, reduc3 )
local coherences = coherence_compute.compute_coherences_self_task(update_neighbours, parts1)

local __demand(__leaf) task self_task([parts1], [config]) where
   [read1_privs], [readconf_privs], [write1_privs], [reduc1_privs], [reducconf_privs], reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}), 
                                  reads(parts1.neighbour_part_space._valid), reads(config.space),
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
               [kernel_name(rexpr [parts1][part1] end, rexpr [parts1][part2] end, rexpr r2 end, rexpr config[0] end)]
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
local read1, read2, read3, write1, write2, write3, reduc1, reduc2, reduc3 = compute_privileges.three_region_privileges( kernel_name )
local cell_pair_task, update_neighbours1 = generate_symmetric_pairwise_task( kernel_name, read1, read2, read3, write1, write2, reduc1, reduc2, reduc3 )
local cell_self_task, update_neighbours2 = generate_symmetric_self_task( kernel_name, read1, read2, read3, write1, write2, reduc1, reduc2, reduc3 )
local update_neighbours = update_neighbours1 or update_neighbours2
local update_neighbours_quote = rquote

end
if update_neighbours then
  local temp_variables = {}
  temp_variables.config = config
  update_neighbours_quote = neighbour_init.update_cells(temp_variables)
end
local start_timing_quote = rquote

end
local end_timing_quote = rquote

end
if DSL_settings.TIMING then
  local starttime = regentlib.newsymbol()
  local endtime = regentlib.newsymbol()
  start_timing_quote = rquote
    var [starttime] = c.legion_get_current_time_in_micros()
  end
  end_timing_quote = rquote
    var [endtime] = c.legion_get_current_time_in_micros();
    [variables.config][0].timing_config.user_task_time = [variables.config][0].timing_config.user_task_time + ([endtime] - [starttime])
  end
end

assert(false, "No High Performance implementation for symmetric kernels exists yet.")

return symmetric

end

function old_create_asymmetric_pairwise_runner( kernel_name, config, cell_space )
local read1, read2, read3, write1, write2, write3, reduc1, reduc2, reduc3 = compute_privileges.three_region_privileges( kernel_name )
--While write2 is computed, asymmetric kernels are not allowed to write to write2
local cell_pair_task, update_neighbours1 = generate_asymmetric_pairwise_task( kernel_name, read1, read2, read3, write1, reduc1, reduc3 )
local cell_self_task, update_neighbours2 = generate_asymmetric_self_task( kernel_name, read1, read2, read3, write1, reduc1, reduc3 )
local update_neighbours = update_neighbours1 or update_neighbours2
local update_neighbours_quote = rquote

end
if update_neighbours then
  local temp_variables = {}
  temp_variables.config = config
  update_neighbours_quote = neighbour_init.update_cells(temp_variables)
end

local start_timing_quote = rquote

end
local end_timing_quote = rquote

end
if DSL_settings.TIMING then
  local starttime = regentlib.newsymbol()
  local endtime = regentlib.newsymbol()
  start_timing_quote = rquote
    var [starttime] = c.legion_get_current_time_in_micros()
  end
  end_timing_quote = rquote
    var [endtime] = c.legion_get_current_time_in_micros();
    [variables.config][0].timing_config.user_task_time = [variables.config][0].timing_config.user_task_time + ([endtime] - [starttime])
  end
end

local asymmetric = rquote
    [start_timing_quote];
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
    --Exit gracefully if max_cutoff >= half the boxsize
    regentlib.assert(x_radii <= x_count / 2, "max_cutoff is too large relative to the boxsize for High Performance implementation. Please use a basic neighbour algorithm instead")
    regentlib.assert(y_radii <= y_count / 2, "max_cutoff is too large relative to the boxsize for High Performance implementation. Please use a basic neighbour algorithm instead")
    regentlib.assert(z_radii <= z_count / 2, "max_cutoff is too large relative to the boxsize for High Performance implementation. Please use a basic neighbour algorithm instead")
--    __demand(__trace)
    do
        for x = -x_radii, x_radii+1 do
            for y = -y_radii, y_radii+1 do
                for z = -z_radii, z_radii + 1 do
                    if x~= 0 or y~=0 or z~= 0 then
                        __demand(__index_launch)
                        for cell in cell_space.colors do
                          --Launch all tasks. Not yet supporting complex large cutoff radii.
                          cell_pair_task(cell_space[cell], cell_space[int3d({wrap(cell.x+x,x_count), wrap(cell.y+y, y_count), wrap(cell.z+z, z_count) })], config)
                        end
                    else
                        __demand(__index_launch)
                        for cell1 in cell_space.colors do
                            cell_self_task(cell_space[cell1], config)
                        end
                    end
                end
            end
        end
    end
    [end_timing_quote];
    [update_neighbours_quote];
end
return asymmetric
end


-----------------------------------------
--End of Launcher Code  -----------------
-----------------------------------------

--Generate a task to be executed on every particle in the system
function generate_per_part_task( kernel_name, read1, write1, reduc1, readconf, reducconf)
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local config = regentlib.newsymbol(region(ispace(int1d), config_type), "config")
--local update_neighbours, read1_privs, write1_privs, reduc1_privs = privilege_lists.get_privileges_self_task( parts1, read1, {}, write1, {}, reduc1, {} )
local update_neighbours, read1_privs, write1_privs, reduc1_privs, readconf_privs, reducconf_privs = 
                                            privilege_lists.get_privileges_per_part(parts1, read1, write1, reduc1, config, readconf, reducconf)
local coherences = coherence_compute.compute_coherences_self_task(update_neighbours, parts1)

local __demand(__leaf) task per_part_task([parts1], [config]) where
   [read1_privs], [write1_privs],[reduc1_privs], [readconf_privs], [reducconf_privs], [coherences], 
   reads( parts1.neighbour_part_space._valid )
   do
   for part1 in [parts1].ispace do
     if [parts1][part1].neighbour_part_space._valid then
         [kernel_name(rexpr [parts1][part1] end, rexpr [config][0] end)]
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
local per_part_task, update_neighbours = generate_per_part_task( kernel_name, read1, write1, reduc1, read2, reduc2 )
local update_neighbours_quote = rquote

end
if update_neighbours then
  local temp_variables = {}
  temp_variables.config = config
  update_neighbours_quote = neighbour_init.update_cells(temp_variables)
end

local start_timing_quote = rquote

end
local end_timing_quote = rquote

end
if DSL_settings.TIMING then
  local starttime = regentlib.newsymbol()
  local endtime = regentlib.newsymbol()
  start_timing_quote = rquote
    var [starttime] = c.legion_get_current_time_in_micros()
  end
  end_timing_quote = rquote
    c.legion_runtime_issue_execution_fence(__runtime(), __context())
    var [endtime] = c.legion_get_current_time_in_micros();
    [variables.config][0].timing_config.user_task_time = [variables.config][0].timing_config.user_task_time + ([endtime] - [starttime])
  end
end
local runner = rquote

    [start_timing_quote];
    --For each cell, call the task!
    __demand(__index_launch)
    for slice in [neighbour_init.supercell_partition].colors do
--    for cell1 in cell_space.colors do
       per_part_task([neighbour_init.supercell_partition][slice], config)
    end
   [end_timing_quote];
   [update_neighbours_quote];
end

return runner
end

-----------------------------------------
--End of Per Part Tasks------------------
-----------------------------------------

--New algorithms for the new neighbour scheme
function new_generate_asymmetric_pairwise_task( kernel_name, read1, read2, read3, write1, reduc1, reduc3 )
--Asymmetric kernel can only write to part1
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local config = regentlib.newsymbol(region(ispace(int1d), config_type), "config")
local allparts = regentlib.newsymbol(region(ispace(int1d), part), "allparts")
local cell_partition = regentlib.newsymbol(partition(disjoint, allparts, ispace(int3d)), "cell_partition")
local  update_neighbours, read1_privs, read2_privs, read3_privs, write1_privs, reduc1_privs, reduc3_privs =
      privilege_lists.get_privileges_asymmetric_pair_task( parts1, parts2, config, read1, read2, read3,
              write1, reduc1, reduc3 )
local coherences = coherence_compute.compute_coherences_pair_task(update_neighbours, parts1, parts2)

local directions = terralib.newlist({
    rexpr int3d({ 0,  1,  0}) end,
    rexpr int3d({ 0,  0,  1}) end,
    rexpr int3d({ 0,  1,  1}) end,
    rexpr int3d({ 0,  1, -1}) end,
    rexpr int3d({ 1,  0,  0}) end,
    rexpr int3d({ 1,  0,  1}) end,
    rexpr int3d({ 1,  1,  0}) end,
    rexpr int3d({ 1,  0, -1}) end,
    rexpr int3d({ 1, -1,  0}) end,
    rexpr int3d({ 1,  1,  1}) end,
    rexpr int3d({ 1, -1, -1}) end,
    rexpr int3d({ 1, -1,  1}) end,
    rexpr int3d({ 1,  1, -1}) end,

    rexpr int3d({ 0, -1,  0}) end,
    rexpr int3d({ 0,  0, -1}) end,
    rexpr int3d({ 0, -1, -1}) end,
    rexpr int3d({ 0, -1,  1}) end,
    rexpr int3d({-1,  0,  0}) end,
    rexpr int3d({-1,  0, -1}) end,
    rexpr int3d({-1, -1,  0}) end,
    rexpr int3d({-1,  0,  1}) end,
    rexpr int3d({-1,  1,  0}) end,
    rexpr int3d({-1, -1, -1}) end,
    rexpr int3d({-1,  1,  1}) end,
    rexpr int3d({-1,  1, -1}) end,
    rexpr int3d({-1, -1,  1}) end
})

local cell = regentlib.newsymbol(int3d, "cell")
local count_xcells = regentlib.newsymbol(int, "count_xcells")
local count_ycells = regentlib.newsymbol(int, "count_ycells")
local count_zcells = regentlib.newsymbol(int, "count_zcells")
local xlo = regentlib.newsymbol()
local xhi = regentlib.newsymbol()
local ylo = regentlib.newsymbol()
local yhi = regentlib.newsymbol()
local zlo = regentlib.newsymbol()
local zhi = regentlib.newsymbol()
local box_x = regentlib.newsymbol()
local box_y = regentlib.newsymbol()
local box_z = regentlib.newsymbol()
local half_box_x = regentlib.newsymbol()
local half_box_y = regentlib.newsymbol()
local half_box_z = regentlib.newsymbol()
local ne_cell = regentlib.newsymbol(int3d)
local nano_start = regentlib.newsymbol(int64)
local nano_end = regentlib.newsymbol(int64)
local nano_total = regentlib.newsymbol(int64)
            --Loop over all directions
--local function interact_loop_func()
--local __quotes = terralib.newlist()
--              for i = 1, 26 do
--                __quotes:insert(rquote
--                    [ne_cell] = ([cell] + [directions[i]] + {[count_xcells],[count_ycells],[count_zcells]})
--                                                           %{[count_xcells],[count_ycells],[count_zcells]}
--                    --Check if the cell is in the supercell or not
--                    if [ne_cell].x < [xlo] or [ne_cell].x >= [xhi] or [ne_cell].y < [ylo] or [ne_cell].y >= [yhi] or [ne_cell].z < [zlo] or [ne_cell].z >= [zhi] then
--                        --Neighbour cell is inside, lets interact the particles!
--                        for part1 in cell_partition[cell].ispace do
--                            if [parts1][part1].neighbour_part_space._valid then
--                                for part2 in cell_partition[ne_cell].ispace do
--                                    if [parts2][part2].neighbour_part_space._valid then
--                                        --Compute the distance between them
--                                        var dx = [parts1][part1].core_part_space.pos_x - [parts2][part2].core_part_space.pos_x
--                                        var dy = [parts1][part1].core_part_space.pos_y - [parts2][part2].core_part_space.pos_y
--                                        var dz = [parts1][part1].core_part_space.pos_z - [parts2][part2].core_part_space.pos_z
--                                        if (dx > [half_box_x]) then dx = dx - [box_x] end
--                                        if (dy > [half_box_y]) then dy = dy - [box_y] end
--                                        if (dz > [half_box_z]) then dz = dz - [box_z] end
--                                        if (dx <-[half_box_x]) then dx = dx + [box_x] end
--                                        if (dy <-[half_box_y]) then dy = dy + [box_y] end
--                                        if (dz <-[half_box_z]) then dz = dz + [box_z] end
--                                        var cutoff2 = regentlib.fmax([parts1][part1].core_part_space.cutoff, [parts2][part2].core_part_space.cutoff)
--                                        cutoff2 = cutoff2 * cutoff2
--                                        var r2 = dx*dx + dy*dy + dz*dz
--                                        if(r2 <= cutoff2) then
--                                          [nano_start] = regentlib.c.legion_get_current_time_in_nanos();
--                                          [kernel_name(rexpr [parts1][part1] end, rexpr [parts2][part2] end, rexpr r2 end, rexpr config[0] end)];
--                                          [nano_end] = regentlib.c.legion_get_current_time_in_nanos();
--                                          [nano_total] = [nano_total] + ([nano_end] - [nano_start])
--                                        end
--                                    end
--                                end
--                            end
--                        end
--                    end
--                end)
--              end
--              return __quotes
--end
--
--local interact_quote = interact_loop_func()

local __demand(__leaf) task asym_pairwise_task([parts1], [parts2], [config], [allparts], 
                                                [cell_partition], cell1 : int3d,
                                       full_sort_list : region(ispace(int1d), sorting_ids), --Use to reference the subcell partition
                                       sort_subcell_partition : partition(disjoint, full_sort_list, ispace(int3d)),
                                       supercell_sort_list : region(ispace(int1d), sorting_ids) ) where 
  [read1_privs], [read2_privs], [read3_privs], [write1_privs], [reduc1_privs], [reduc3_privs], reads(config.space, config.neighbour_config),
  reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}),
  reads(parts2.core_part_space.{pos_x, pos_y, pos_z, cutoff}), reads(parts1.neighbour_part_space._valid), reads(parts2.neighbour_part_space._valid),
  reads(parts1.neighbour_part_space.sorting_positions), reads(parts2.neighbour_part_space.sorting_positions),
  [coherences], reads(supercell_sort_list), reads(full_sort_list), reads(allparts.neighbour_part_space),
--TODO REMOVE
   reads(parts2.neighbour_part_space.cell_id) do

    --parts1 is our supercell we write to. parts2 is our halo we're reading from. In this task we only interact subcells in our supercell with subcells from the halo
    --This could change in the future, but for now makes life easier
    var x_per_super = config[0].neighbour_config.x_cells / config[0].neighbour_config.x_supercells
    var y_per_super = config[0].neighbour_config.y_cells / config[0].neighbour_config.y_supercells
    var z_per_super = config[0].neighbour_config.z_cells / config[0].neighbour_config.z_supercells
    var [xlo] = (cell1.x) * x_per_super
    var [xhi] = (cell1.x+1) * x_per_super --Loops non inclusive
    var [ylo] = cell1.y * y_per_super
    var [yhi] = (cell1.y+1) * y_per_super
    var [zlo] = (cell1.z) * z_per_super
    var [zhi] = (cell1.z+1) * z_per_super

    var [count_xcells] = config[0].neighbour_config.x_cells
    var [count_ycells] = config[0].neighbour_config.y_cells
    var [count_zcells] = config[0].neighbour_config.z_cells
    var [box_x] = config[0].space.dim_x
    var [box_y] = config[0].space.dim_y
    var [box_z] = config[0].space.dim_z
    var [half_box_x] = 0.5 * box_x
    var [half_box_y] = 0.5 * box_y
    var [half_box_z] = 0.5 * box_z

    regentlib.assert(config[0].neighbour_config.cell_dim_x > config[0].neighbour_config.max_cutoff, 
                    "Cells couldn't be created small enough to support High Performance implementation")
    regentlib.assert(config[0].neighbour_config.cell_dim_y > config[0].neighbour_config.max_cutoff, 
                    "Cells couldn't be created small enough to support High Performance implementation")
    regentlib.assert(config[0].neighbour_config.cell_dim_z > config[0].neighbour_config.max_cutoff, 
                    "Cells couldn't be created small enough to support High Performance implementation")

    var starttime = c.legion_get_current_time_in_micros()
    var [ne_cell];
    var [nano_start];
    var [nano_end];
    var [nano_total] = 0;
    var nano_total2 : int64 = 0;

    var total = 0;
    var hits = 0;

    var direction_array : int3d[26]
    var periodic_sort_distance : double[26]

    [(function() local __quotes = terralib.newlist()
        for i=1, 26 do
            __quotes:insert(rquote
                direction_array[i-1] = [directions[i]]
            end)
        end                                                                                                                                                                        return __quotes
    end) ()];
    for i=0, 26 do
        var periodic_sort_x = double(direction_array[i].x) * [box_x]
        var periodic_sort_y = double(direction_array[i].y) * [box_y]
        var periodic_sort_z = double(direction_array[i].z) * [box_z]
        periodic_sort_distance[i] = periodic_sort_x * periodic_sort_x + periodic_sort_y * periodic_sort_y + periodic_sort_z * periodic_sort_z
        periodic_sort_distance[i] = sqrt(periodic_sort_distance[i]) 
    end

    var max_cutoff = config[0].neighbour_config.max_cutoff

    --For now implement this in a bad way
    for x=xlo, xhi do
    for y=ylo, yhi do
    for z=zlo, zhi do

        var [cell] = int3d({x, y, z})
        --Find my bounds
        var cell_hi = sort_subcell_partition[cell].ispace.bounds.hi
        var cell_lo = sort_subcell_partition[cell].ispace.bounds.lo
        for i = int(cell_lo), int(cell_hi)+1 do
            if supercell_sort_list[int1d(i)].sid[0] == int1d(-1) then
                cell_hi = int1d(i)
                break
            end
        end 
        --Do positive directions first
        for dir = 0, 13 do
            var ne_cell : int3d = (cell + direction_array[dir] +{count_xcells,count_ycells,count_zcells})%{count_xcells,count_ycells,count_zcells}
            --Check if the cell is in the supercell or not
            if [ne_cell].x < [xlo] or [ne_cell].x >= [xhi] or [ne_cell].y < [ylo] or [ne_cell].y >= [yhi] or [ne_cell].z < [zlo] or [ne_cell].z >= [zhi] then
                --Neighbour cell is a halo, lets interact the particles!
                
                --Find neighbour cell bounds
                var necell_hi = sort_subcell_partition[ne_cell].ispace.bounds.hi + 1
                var necell_lo = sort_subcell_partition[ne_cell].ispace.bounds.lo
                --Loop over our particles in reverse order
                for i= int(cell_hi-1), int(cell_lo)-1, -1 do
                    --Find our particle index
                    var part1 = supercell_sort_list[int1d(i)].sid[dir]
                    --Loop over our neighbour in ascending order
                    for j = int(necell_lo), int(necell_hi) do
                        var part2 = full_sort_list[int1d(j)].sid[dir]
                        --Check if its a valid index
                        if part2 == int1d(-1) then
                            --If not we reset necell_hi and stop
                            necell_hi = int1d(j)
                            break
                        end
                        --Its a valid one, check for the sorted distance
                        var sort_distance = [parts1][part1].neighbour_part_space.sorting_positions[dir] 
                                          - [parts2][part2].neighbour_part_space.sorting_positions[dir]
                        if abs(sort_distance) > 0.5 * periodic_sort_distance[dir] then
                            sort_distance = abs(sort_distance) - periodic_sort_distance[dir]
                        end
--                        if [parts1][part1].core_part_space.id == int1d(3) and ne_cell == int3d({0, 2, 1}) then
--                            format.println("Checking 3 and {}, cell {} {} {}", [parts2][part2].core_part_space.id, ne_cell.x, ne_cell.y, ne_cell.z)
--                            for k = int(necell_lo), int(necell_hi) do
--                                format.println("Not yet checked {}", [parts2][full_sort_list[int1d(k)].sid[dir]].core_part_space.id)
--                                var temp_sd = [parts1][part1].neighbour_part_space.sorting_positions[dir]
--                                            - [parts2][full_sort_list[int1d(k)].sid[dir]].neighbour_part_space.sorting_positions[dir]
--                                if abs(temp_sd) >  0.5 * periodic_sort_distance[dir] then
--                                    temp_sd = abs(temp_sd) - periodic_sort_distance[dir]
--                                end
--                                format.println("Sort distance would be {} and cutoff is {}", temp_sd, max_cutoff)
--                            end
--                        end
                        --If its out of range, reset necell_hi and stop
                        if abs(sort_distance) > max_cutoff then
                            necell_hi = int1d(j)
                            break
                        end
                        --This pair can be in range, so do the real computation
                        --Compute the distance between them
                        total = total + 1;
                        nano_start = regentlib.c.legion_get_current_time_in_nanos();
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
                        nano_end = regentlib.c.legion_get_current_time_in_nanos();
                        nano_total2 = nano_total2 + (nano_end - nano_start)
                        if(r2 <= cutoff2) then
                          nano_start = regentlib.c.legion_get_current_time_in_nanos();
                          [kernel_name(rexpr [parts1][part1] end, rexpr [parts2][part2] end, rexpr r2 end, rexpr config[0] end)];
                          nano_end = regentlib.c.legion_get_current_time_in_nanos();
                          nano_total = nano_total + (nano_end - nano_start)
                          hits = hits + 1
                        end
                    end
                    
                end
                
            end
        end
        for dir = 13, 26 do
            var ne_cell : int3d = (cell + direction_array[dir] +{count_xcells,count_ycells,count_zcells})%{count_xcells,count_ycells,count_zcells}
            --Check if the cell is in the supercell or not
            if [ne_cell].x < [xlo] or [ne_cell].x >= [xhi] or [ne_cell].y < [ylo] or [ne_cell].y >= [yhi] or [ne_cell].z < [zlo] or [ne_cell].z >= [zhi] then
                --Neighbour cell is a halo, lets interact the particles!
                
                --Find neighbour cell bounds
                var necell_hi = sort_subcell_partition[ne_cell].ispace.bounds.hi
                var necell_lo = sort_subcell_partition[ne_cell].ispace.bounds.lo
                for i = int(cell_lo), int(cell_hi) do
                    var part1 = supercell_sort_list[int1d(i)].sid[dir-13]
                    for j = int(necell_hi), int(necell_lo)-1, -1 do
                       var part2 = full_sort_list[int1d(j)].sid[dir-13] 
                       --Check if its a valid index
                       if part2 == int1d(-1) then
                           --Its not a valid one
                           necell_hi = necell_hi - 1
                       else
                           --Its a valid one, check for the sorted distance
                           var sort_distance = [parts1][part1].neighbour_part_space.sorting_positions[dir] 
                                             - [parts2][part2].neighbour_part_space.sorting_positions[dir]
                            if abs(sort_distance) > 0.5 * periodic_sort_distance[dir] then
                                sort_distance = abs(sort_distance) - periodic_sort_distance[dir]
                            end
                          if [parts1][part1].core_part_space.id == int1d(3) and ne_cell == int3d({2, 2, 1}) then
                                 format.println("Checking 3 and {}, cell {} {} {}", [parts2][part2].core_part_space.id, ne_cell.x, ne_cell.y, ne_cell.z)
                            for k = int(necell_hi), int(necell_lo)-1, -1 do
                                format.println("Not yet checked {}", [parts2][full_sort_list[int1d(k)].sid[dir]].core_part_space.id)
                                if [parts2][full_sort_list[int1d(k)].sid[dir]].core_part_space.id == int1d(994) then
                                    var temp_sd = [parts1][part1].neighbour_part_space.sorting_positions[dir]
                                                + [parts2][full_sort_list[int1d(k)].sid[dir]].neighbour_part_space.sorting_positions[dir]
                                    format.println("Direction is {}", dir)
                                    format.println("3 pos is {} {} {}, sort pos is {}", [parts1][part1].core_part_space.pos_x, 
                                                    [parts1][part1].core_part_space.pos_y, [parts1][part1].core_part_space.pos_z,
                                                    [parts1][part1].neighbour_part_space.sorting_positions[dir])
                                    format.println("94 pos is {} {} {}, sort pos is {}", [parts2][full_sort_list[int1d(k)].sid[dir]].core_part_space.pos_x,
                                     [parts2][full_sort_list[int1d(k)].sid[dir]].core_part_space.pos_y,
                                     [parts2][full_sort_list[int1d(k)].sid[dir]].core_part_space.pos_z,  [parts2][full_sort_list[int1d(k)].sid[dir]].neighbour_part_space.sorting_positions[dir])
                                    format.println("vector is {} {} {}", direction_array[dir].x, direction_array[dir].y, direction_array[dir].z)
                                    format.println("Non-wrapped sort distance is {}", temp_sd)
                                    if abs(temp_sd) >  0.5 * periodic_sort_distance[dir] then
                                        temp_sd = abs(temp_sd) - periodic_sort_distance[dir]
                                    end
                                    format.println("Sort distance would be {} and cutoff is {}", temp_sd, max_cutoff)
                                    var dx = [parts1][part1].core_part_space.pos_x - [parts2][full_sort_list[int1d(k)].sid[dir]].core_part_space.pos_x
                                    var dy = [parts1][part1].core_part_space.pos_y - [parts2][full_sort_list[int1d(k)].sid[dir]].core_part_space.pos_y
                                    var dz = [parts1][part1].core_part_space.pos_z - [parts2][full_sort_list[int1d(k)].sid[dir]].core_part_space.pos_z
                                    if (dx > half_box_x) then dx = dx - box_x end
                                    if (dy > half_box_y) then dy = dy - box_y end
                                    if (dz > half_box_z) then dz = dz - box_z end
                                    if (dx <-half_box_x) then dx = dx + box_x end
                                    if (dy <-half_box_y) then dy = dy + box_y end
                                    if (dz <-half_box_z) then dz = dz + box_z end
                                    var cutoff2 = regentlib.fmax([parts1][part1].core_part_space.cutoff, [parts2][full_sort_list[int1d(k)].sid[dir]].core_part_space.cutoff)
                                    cutoff2 = cutoff2 * cutoff2
                                    var r2 = dx*dx + dy*dy + dz*dz
                                    format.println("Part 3 @ {} {} {}", [parts1][part1].core_part_space.pos_x, [parts1][part1].core_part_space.pos_y, 
                                                                        [parts1][part1].core_part_space.pos_z)
                                    format.println("Part 94 @ {} {} {}", [parts2][full_sort_list[int1d(k)].sid[dir]].core_part_space.pos_x, [parts1][full_sort_list[int1d(k)].sid[dir]].core_part_space.pos_y, 
                                                                        [parts2][full_sort_list[int1d(k)].sid[dir]].core_part_space.pos_z)

                                    format.println("{} {} {}, Real distance2 {}, cutoff2 = {}",dx, dy, dz, r2, cutoff2)
                                end
                            end
                                 format.println("Sort distance {}, max_cutoff {}", sort_distance, max_cutoff)
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
                          end
                           --If its out of range, reset necell_lo and stop
                           if abs(sort_distance) > max_cutoff then
                               necell_lo = int1d(j)
                               break
                           end
                           --This pair can be in range, so do the real computation
                           --Compute the distance between them
                           total = total + 1;
                           nano_start = regentlib.c.legion_get_current_time_in_nanos();
                           regentlib.assert(int1d(part1) ~= int1d(part2), "part1 == part2") --TODO
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
                           nano_end = regentlib.c.legion_get_current_time_in_nanos();
                           nano_total2 = nano_total2 + (nano_end - nano_start)
                           if(r2 <= cutoff2) then
                             nano_start = regentlib.c.legion_get_current_time_in_nanos();
                             [kernel_name(rexpr [parts1][part1] end, rexpr [parts2][part2] end, rexpr r2 end, rexpr config[0] end)];
                             nano_end = regentlib.c.legion_get_current_time_in_nanos();
                             nano_total = nano_total + (nano_end - nano_start)
                             hits = hits + 1
                           end
                       end
                    end
                end
            end
        end

    end
    end
    end

--    for z = zlo, zhi do
--        --First loop, loop over cells at xlo and xhi -- this covers the "right" and and "left" surfaces of the cell
--        for y = ylo, yhi do
--            --xlo cell
--            var [cell] = int3d({xlo, y, z});
--            [interact_quote];
--
--            --xhi -1 cell
--            [cell] = int3d({xhi-1, y, z});
--            --Loop over all directions
--            [interact_quote];
--            
--        end --loop from ylo to yhi
--       
--        --Next we do xlo+1 to xhi-1 for the ylo and yhi cells -- this covers the "top" and "bottom" surfaces of the cell that weren't covered by the x cells
--        for x = xlo+1, xhi-1 do
--            --ylo cell
--            var [cell] = int3d({x, ylo, z});
--            --Loop over all directions
--            [interact_quote];
--            
--            --yhi-1 cell
--            [cell] = int3d({x, yhi-1, z});
--            --Loop over all directions
--            [interact_quote];
--        end
--    end --end of zcell
--
--    --Finally loop over the front and back surfaces not covered by the other loops
--    for x= xlo+1, xhi-1 do
--        for y=ylo+1, yhi-1 do
--            --zlo cell
--            var [cell] = int3d({x, y, zlo});
--            [interact_quote]
--
--            --zhi-1 cell
--            cell = int3d({x, y, zhi-1});
--            [interact_quote]
--        end
--    end

    --Should have computed all interactions with the halos now
    var endtime = c.legion_get_current_time_in_micros()

--    format.println("PAIR: runtime was {}us, interaction time was {}us, distance time was {}us", (endtime-starttime), (nano_total / 1000), (nano_total2 / 1000))
--    format.println("PAIR: hits {} misses {}", hits, total-hits)
    --format.println("runtime was {}us, interaction time was {}us", (endtime-starttime), ([nano_total] / 1000))
end
return asym_pairwise_task, update_neighbours
end



--Generate a self task
--This function assumes the cutoff of only the updated part is relevant
function new_generate_asymmetric_self_task( kernel_name, read1, read2, read3, write1, reduc1, reduc3 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local config = regentlib.newsymbol(region(ispace(int1d), config_type), "config")
local update_neighbours, read1_privs, readconf_privs, write1_privs, reduc1_privs, reducconf_privs = 
                                                                   privilege_lists.get_privileges_self_task( parts1, config, read1, read2, read3,
                                                                                                             write1, {}, reduc1, {}, reduc3 )
local coherences = coherence_compute.compute_coherences_self_task(update_neighbours, parts1)

local directions = terralib.newlist({
    rexpr int3d({ 0,  1,  0}) end,
    rexpr int3d({ 0,  0,  1}) end,
    rexpr int3d({ 0,  1,  1}) end,
    rexpr int3d({ 0,  1, -1}) end,
    rexpr int3d({ 1,  0,  0}) end,
    rexpr int3d({ 1,  0,  1}) end,
    rexpr int3d({ 1,  1,  0}) end,
    rexpr int3d({ 1,  0, -1}) end,
    rexpr int3d({ 1, -1,  0}) end,
    rexpr int3d({ 1,  1,  1}) end,
    rexpr int3d({ 1, -1, -1}) end,
    rexpr int3d({ 1, -1,  1}) end,
    rexpr int3d({ 1,  1, -1}) end,

    rexpr int3d({ 0, -1,  0}) end,
    rexpr int3d({ 0,  0, -1}) end,
    rexpr int3d({ 0, -1, -1}) end,
    rexpr int3d({ 0, -1,  1}) end,
    rexpr int3d({-1,  0,  0}) end,
    rexpr int3d({-1,  0, -1}) end,
    rexpr int3d({-1, -1,  0}) end,
    rexpr int3d({-1,  0,  1}) end,
    rexpr int3d({-1,  1,  0}) end,
    rexpr int3d({-1, -1, -1}) end,
    rexpr int3d({-1,  1,  1}) end,
    rexpr int3d({-1,  1, -1}) end,
    rexpr int3d({-1, -1,  1}) end
})

local __demand(__leaf) task self_task([parts1], [config],allparts : region(ispace(int1d), part),
                                       cell_partition : partition(disjoint, allparts, ispace(int3d)), cell1 : int3d,
                                       full_sort_list : region(ispace(int1d), sorting_ids), --Use to reference the subcell partition
                                       sort_subcell_partition : partition(disjoint, full_sort_list, ispace(int3d)),
                                       supercell_sort_list : region(ispace(int1d), sorting_ids) ) where
   [read1_privs], [readconf_privs], [write1_privs], [reduc1_privs], [reducconf_privs], reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}), 
                                  reads(parts1.neighbour_part_space._valid, parts1.neighbour_part_space.sorting_positions), 
                                  reads(config.space, config.neighbour_config),
                                  reads(supercell_sort_list),
   [coherences] do

    var x_per_super = config[0].neighbour_config.x_cells / config[0].neighbour_config.x_supercells
    var y_per_super = config[0].neighbour_config.y_cells / config[0].neighbour_config.y_supercells
    var z_per_super = config[0].neighbour_config.z_cells / config[0].neighbour_config.z_supercells
    var xlo = (cell1.x) * x_per_super
    var xhi = (cell1.x+1) * x_per_super --Loops non inclusive
    var ylo = cell1.y * y_per_super
    var yhi = (cell1.y+1) * y_per_super
    var zlo = (cell1.z) * z_per_super
    var zhi = (cell1.z+1) * z_per_super

    var count_xcells = config[0].neighbour_config.x_cells
    var count_ycells = config[0].neighbour_config.y_cells
    var count_zcells = config[0].neighbour_config.z_cells
    var box_x = config[0].space.dim_x
    var box_y = config[0].space.dim_y
    var box_z = config[0].space.dim_z
    var half_box_x = 0.5 * box_x
    var half_box_y = 0.5 * box_y
    var half_box_z = 0.5 * box_z

    var direction_array : int3d[13]

    [(function() local __quotes = terralib.newlist()
        for i=1, 13 do
            __quotes:insert(rquote
                direction_array[i-1] = [directions[i]]
            end)
        end                                                                                                                                                                        return __quotes
    end) ()];

    regentlib.assert(config[0].neighbour_config.cell_dim_x > config[0].neighbour_config.max_cutoff, 
                    "Cells couldn't be created small enough to support High Performance implementation")
    regentlib.assert(config[0].neighbour_config.cell_dim_y > config[0].neighbour_config.max_cutoff, 
                    "Cells couldn't be created small enough to support High Performance implementation")
    regentlib.assert(config[0].neighbour_config.cell_dim_z > config[0].neighbour_config.max_cutoff, 
                    "Cells couldn't be created small enough to support High Performance implementation")

    var starttime = c.legion_get_current_time_in_micros()
    var nano_start : int64;
    var nano_end : int64;
    var nano_total : int64 = 0;
    var nano_total2 : int64 = 0;
    var total = 0;
    var hits = 0;
    var max_cutoff = config[0].neighbour_config.max_cutoff
    --Loop over all the internal cells
    for x = xlo, xhi do
        for y = ylo, yhi do
            for z = zlo, zhi do
                --get the cell
                var cell : int3d = int3d({x, y, z});
                --Cell with itself (self task section)
                for part1 in cell_partition[cell].ispace do
                    if [parts1][part1].neighbour_part_space._valid then
                        for part2 in cell_partition[cell].ispace do
                            if int1d(part1) ~= int1d(part2) and [parts1][part2].neighbour_part_space._valid then
                                --Compute the distance between them
                                total = total + 1;
                                nano_start = regentlib.c.legion_get_current_time_in_nanos();
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
                                nano_end = regentlib.c.legion_get_current_time_in_nanos();
                                nano_total2 = nano_total2 + (nano_end - nano_start)
--                            if [parts1][part1].core_part_space.id == int1d(0) then
--                                 format.println("Checking 0 and {}, cell {} {} {}", [parts1][part2].core_part_space.id, cell.x, cell.y, cell.z)
--                            end
                                if(r2 <= cutoff2) then
                                  nano_start = regentlib.c.legion_get_current_time_in_nanos();
                                  [kernel_name(rexpr [parts1][part1] end, rexpr [parts1][part2] end, rexpr r2 end, rexpr config[0] end)];
                                  nano_end = regentlib.c.legion_get_current_time_in_nanos();
                                  nano_total = nano_total + (nano_end - nano_start)
                                  hits = hits + 1
                                end
                            end
                        end
                    end
                end
                --Find my bounds
                var cell_hi = sort_subcell_partition[cell].ispace.bounds.hi
                var cell_lo = sort_subcell_partition[cell].ispace.bounds.lo
                for i = int(cell_lo), int(cell_hi)+1 do
                    if supercell_sort_list[int1d(i)].sid[0] == int1d(-1) then
                        cell_hi = int1d(i)
                        break
                    end
                end
                cell_hi = cell_hi - 1 
                --Loop over all neighbouring cells in the positive direction, if in range we do both directions
                for dir = 0, 13 do
                    var ne_cell : int3d = (cell + direction_array[dir] +{count_xcells,count_ycells,count_zcells})%{count_xcells,count_ycells,count_zcells}
                    --Check if the cell is in the supercell or not
                    if ne_cell.x >= xlo and ne_cell.x < xhi and ne_cell.y >= ylo and ne_cell.y < yhi and ne_cell.z >= zlo and ne_cell.z < zhi then
                        --Neighbour cell is inside, lets interact the particles!
                        
                        --Find neighbour cell bounds
                        var necell_hi = sort_subcell_partition[ne_cell].ispace.bounds.hi + 1
                        var necell_lo = sort_subcell_partition[ne_cell].ispace.bounds.lo
    
                        --Loop over our particles in reverse order
                        for i= int(cell_hi), int(cell_lo)-1, -1 do
                            --Find our particle index
                            var part1 = supercell_sort_list[int1d(i)].sid[dir]
                            --Loop over our neighbour in ascending order
                            for j = int(necell_lo), int(necell_hi) do
                                var part2 = supercell_sort_list[int1d(j)].sid[dir]
                                --Check if its a valid index
                                if part2 == int1d(-1) then
                                    --If not we reset necell_hi and stop
                                    necell_hi = int1d(j)
                                    break
                                end
                                --Its a valid one, check for the sorted distance
                                var sort_distance = [parts1][part1].neighbour_part_space.sorting_positions[dir] 
                                                  - [parts1][part2].neighbour_part_space.sorting_positions[dir]
--                            if [parts1][part1].core_part_space.id == int1d(0) then
--                                 format.println("Checking 0 and {}, cell {} {} {}", [parts1][part2].core_part_space.id, ne_cell.x, ne_cell.y, ne_cell.z)
--                            end
                                --If its out of range, reset necell_hi and stop
                                if abs(sort_distance) > max_cutoff then
                                    --Debug statement
--                                    var dx = [parts1][part1].core_part_space.pos_x - [parts1][part2].core_part_space.pos_x
--                                    var dy = [parts1][part1].core_part_space.pos_y - [parts1][part2].core_part_space.pos_y
--                                    var dz = [parts1][part1].core_part_space.pos_z - [parts1][part2].core_part_space.pos_z
--                                    var cutoff2 = regentlib.fmax([parts1][part1].core_part_space.cutoff, [parts1][part2].core_part_space.cutoff)
--                                    cutoff2 = cutoff2 * cutoff2
--                                    var r2 = dx*dx + dy*dy + dz*dz
--                                    if r2 <= cutoff2 then
--                                        format.println("Sort positions {} {} give distance {}",
--                                                        [parts1][part1].neighbour_part_space.sorting_positions[dir],
--                                                        [parts1][part2].neighbour_part_space.sorting_positions[dir],
--                                                        abs(sort_distance))
--                                        format.println("Real positions ({} {} {}) and ({} {} {})", 
--                                                        [parts1][part1].core_part_space.pos_x,
--                                                        [parts1][part1].core_part_space.pos_y,
--                                                        [parts1][part1].core_part_space.pos_z,
--                                                        [parts1][part1].core_part_space.pos_x,
--                                                        [parts1][part2].core_part_space.pos_y,
--                                                        [parts1][part2].core_part_space.pos_z)
--                                        format.println("Real distance was {}", sqrt(r2))
--                                    end
--                                    regentlib.assert(cutoff2 < r2, "Failed to correctly get neighbours")
                                    --
                                    necell_hi = int1d(j)
                                    break
                                end
                                --This pair can be in range, so do the real computation
                                --Compute the distance between them
                                total = total + 1;
                                nano_start = regentlib.c.legion_get_current_time_in_nanos();
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
                                nano_end = regentlib.c.legion_get_current_time_in_nanos();
                                nano_total2 = nano_total2 + (nano_end - nano_start)
                                if(r2 <= cutoff2) then
                                  nano_start = regentlib.c.legion_get_current_time_in_nanos();
                                  [kernel_name(rexpr [parts1][part1] end, rexpr [parts1][part2] end, rexpr r2 end, rexpr config[0] end)];
                                  [kernel_name(rexpr [parts1][part2] end, rexpr [parts1][part1] end, rexpr r2 end, rexpr config[0] end)];
                                  nano_end = regentlib.c.legion_get_current_time_in_nanos();
                                  nano_total = nano_total + (nano_end - nano_start)
                                  hits = hits + 1
                                end
                            end
                            
                        end
                        
                    end
                end
                --Loop over all neighbouring cells, and if inside the supercell then compute interactions
--                [(function() local __quotes = terralib.newlist()
--                  local ne_cell = regentlib.newsymbol(int3d)
--                  __quotes:insert(rquote
--                        var [ne_cell];
--                    end)
--                  for i = 1, 26 do
--                    __quotes:insert(rquote
--                        [ne_cell] = (cell + [directions[i]] + {count_xcells,count_ycells,count_zcells})%{count_xcells,count_ycells,count_zcells}
--                        --Check if the cell is in the supercell or not
--                        if [ne_cell].x >= xlo and [ne_cell].x < xhi and [ne_cell].y >= ylo and [ne_cell].y < yhi and [ne_cell].z >= zlo and [ne_cell].z < zhi then
--                            --Neighbour cell is inside, lets interact the particles!
--                            for part1 in cell_partition[cell].ispace do
--                                if [parts1][part1].neighbour_part_space._valid then
--                                    for part2 in cell_partition[ne_cell].ispace do
--                                        if [parts1][part2].neighbour_part_space._valid then
--                                            --Compute the distance between them
--                                            total = total + 1;
--                                            nano_start = regentlib.c.legion_get_current_time_in_nanos();
--                                            var dx = [parts1][part1].core_part_space.pos_x - [parts1][part2].core_part_space.pos_x
--                                            var dy = [parts1][part1].core_part_space.pos_y - [parts1][part2].core_part_space.pos_y
--                                            var dz = [parts1][part1].core_part_space.pos_z - [parts1][part2].core_part_space.pos_z
--                                            if (dx > half_box_x) then dx = dx - box_x end
--                                            if (dy > half_box_y) then dy = dy - box_y end
--                                            if (dz > half_box_z) then dz = dz - box_z end
--                                            if (dx <-half_box_x) then dx = dx + box_x end
--                                            if (dy <-half_box_y) then dy = dy + box_y end
--                                            if (dz <-half_box_z) then dz = dz + box_z end
--                                            var cutoff2 = regentlib.fmax([parts1][part1].core_part_space.cutoff, [parts1][part2].core_part_space.cutoff)
--                                            cutoff2 = cutoff2 * cutoff2
--                                            var r2 = dx*dx + dy*dy + dz*dz
--                                            nano_end = regentlib.c.legion_get_current_time_in_nanos();
--                                            nano_total2 = nano_total2 + (nano_end - nano_start)
--                                            if(r2 <= cutoff2) then
--                                              nano_start = regentlib.c.legion_get_current_time_in_nanos();
--                                              [kernel_name(rexpr [parts1][part1] end, rexpr [parts1][part2] end, rexpr r2 end, rexpr config[0] end)];
--                                              nano_end = regentlib.c.legion_get_current_time_in_nanos();
--                                              nano_total = nano_total + (nano_end - nano_start)
--                                              hits = hits + 1
--                                            end
--                                        end
--                                    end
--                                end
--                            end
--                        end
--                    end)
--                  end
--                  return __quotes
--                 end) ()];
            end
        end
    end
    var cell_count = x_per_super * y_per_super * z_per_super
    var endtime = c.legion_get_current_time_in_micros()
--    format.println("SELF: runtime was {}us, interaction time was {}us, distance time was {}us", (endtime-starttime), (nano_total / 1000), (nano_total2 / 1000))
--    format.println("SELF: hits {} misses {}", hits, total-hits)
end
return self_task, update_neighbours
end


--Just need to generate launcher now
function create_asymmetric_pairwise_runner( kernel_name, config, cell_space )
local read1, read2, read3, write1, write2, write3, reduc1, reduc2, reduc3 = compute_privileges.three_region_privileges( kernel_name )
--While write2 is computed, asymmetric kernels are not allowed to write to write2
local cell_halo_task, update_neighbours1 = new_generate_asymmetric_pairwise_task( kernel_name, read1, read2, read3, write1, reduc1, reduc3 )
local cell_self_task, update_neighbours2 = new_generate_asymmetric_self_task( kernel_name, read1, read2, read3, write1, reduc1, reduc3 )
local update_neighbours = update_neighbours1 or update_neighbours2
local update_neighbours_quote = rquote

end
if update_neighbours then
  local temp_variables = {}
  temp_variables.config = config
  update_neighbours_quote = neighbour_init.update_cells(temp_variables)
end

local start_timing_quote = rquote

end
local end_timing_quote = rquote

end
if DSL_settings.TIMING then
  local starttime = regentlib.newsymbol()
  local endtime = regentlib.newsymbol()
  start_timing_quote = rquote
    var [starttime] = c.legion_get_current_time_in_micros()
  end
  end_timing_quote = rquote
    var [endtime] = c.legion_get_current_time_in_micros();
    [variables.config][0].timing_config.user_task_time = [variables.config][0].timing_config.user_task_time + ([endtime] - [starttime])
  end
end

local asymmetric = rquote
    [start_timing_quote];

    --Launch all of the cell + halo tasks
    __demand(__index_launch)
    for supercell in neighbour_init.supercell_partition.colors do
        cell_halo_task(neighbour_init.supercell_partition[supercell], neighbour_init.halo_partition[supercell], variables.config,
                        neighbour_init.padded_particle_array, neighbour_init.cell_partition, supercell,
                        neighbour_init.sorting_id_array,
                        neighbour_init.sorting_array_cell_partition, neighbour_init.sorting_array_supercell_partition[supercell]);
    end

    --Launch all of the cell internal tasks
    __demand(__index_launch)
    for supercell in neighbour_init.supercell_partition.colors do
        cell_self_task( neighbour_init.supercell_partition[supercell], variables.config, neighbour_init.padded_particle_array,
                        neighbour_init.cell_partition, supercell, neighbour_init.sorting_id_array,
                        neighbour_init.sorting_array_cell_partition, neighbour_init.sorting_array_supercell_partition[supercell]);
    end

    [end_timing_quote];
    [update_neighbours_quote];
end
return asymmetric
end

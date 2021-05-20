import "regent"

local compute_privileges = require("src/utils/compute_privilege")
local format = require("std/format")
local string_to_field_path = require("src/utils/string_to_fieldpath")
local coherence_compute = require("src/utils/coherence_computing")
local privilege_lists = require("src/utils/privilege_lists")
local c = regentlib.c



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
   [read1_privs], [write1_privs],[reduc1_privs], [readconf_privs], [reducconf_privs], [coherences]
   do
   for part1 in [parts1].ispace do
      [kernel_name(rexpr [parts1][part1] end, rexpr [config][0] end)]
   end
end
return per_part_task, update_neighbours
end


function run_per_particle_task( kernel_name, config, cell_space )
local read1, read2, write1, write2, reduc1, reduc2 = compute_privileges.two_region_privileges(kernel_name)
local need_tradequeues = false
local per_part_task, update_neighbours = generate_per_part_task( kernel_name, read1, write1, reduc1, read2, reduc2 )

local runner = rquote

   per_part_task(variables.particle_array, config)

end
    return runner
end

function generate_symmetric_self_task( kernel_name, read1, read2, read3, write1, write2, reduc1, reduc2, reduc3 )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local config = regentlib.newsymbol(region(ispace(int1d), config_type), "config")
local update_neighbours, read1_privs, readconf_privs, write1_privs, reduc1_privs, reducconf_privs = privilege_lists.get_privileges_self_task( parts1, config, read1, read2, read3, write1, write2, reduc1, reduc2, reduc3 )
local coherences = coherence_compute.compute_coherences_self_task(update_neighbours, parts1)

local __demand(__leaf) task self_task([parts1], [config]) where
  [read1_privs], [readconf_privs], [write1_privs], [reduc1_privs], [reducconf_privs], reads(parts1.core_part_space.{pos_x, pos_y, pos_z, cutoff}),
                                 reads(config.space),
   [coherences] do
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
               [kernel_name(rexpr [parts1][part1] end, rexpr [parts1][part2] end, rexpr r2 end, rexpr config[0] end)]
             end
           end
       end
   end
end
return self_task, update_neighbours
end


function create_symmetric_pairwise_runner( kernel_name, config, cell_space )
local read1, read2, read3, write1, write2, write3, reduc1, reduc2, reduc3 = compute_privileges.three_region_privileges( kernel_name )
local cell_self_task, update_neighbours = generate_symmetric_self_task( kernel_name, read1, read2, read3, write1, write2, reduc1, reduc2, reduc3 )

local symmetric = rquote
    cell_self_task(variables.particle_array, config)
end
    return symmetric
end

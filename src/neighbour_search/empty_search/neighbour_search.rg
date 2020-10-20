import "regent"

require("defaults")

local compute_privileges = require("src/utils/compute_privilege")
local format = require("std/format")
local string_to_field_path = require("src/utils/string_to_fieldpath")


--Generate the classic MD-style symmetric update pairwise task.
--This function assumes the cutoff is the same for both particles
function generate_symmetric_pairwise_task( kernel_name )
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
  reads(parts2.core_part_space.{pos_x, pos_y, pos_z, cutoff}) do
  for part1 in parts1 do
    for part2 in parts2 do
         [kernel_name(rexpr parts1[part1] end, rexpr parts2[part2] end) ]
    end
  end
end
return pairwise_task, update_neighbours
end

function generate_asymmetric_pairwise_task( kernel_name )
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
  reads(parts2.core_part_space.{pos_x, pos_y, pos_z, cutoff}) do
   for part1 in parts1 do
     for part2 in parts2 do
       [kernel_name( rexpr parts1[part1]end, rexpr parts2[part2] end)]
     end
   end
end
return pairwise_task, update_neighbours
end

--No example self task implemented here

function run_symmetric_pairwise_task( kernel_name )
local read1, read2, write1, write2 = compute_privileges.two_region_privileges( kernel_name )
local cell_pair_task, update_neighbours1 = generate_symmetric_pairwise_task( kernel_name, read1, read2, write1, write2 )
local update_neighbours = update_neighbours1 
local update_neighbours_quote = rquote

end

local symmetric = rquote
   --Tasks to run goes here
  [update_neighbours_quote];
end
return symmetric
end




-----------------------------------------
--End of Pair Tasks ---------------------
-----------------------------------------

--Generate a task to be executed on every particle in the system
function generate_per_part_task( kernel_name )
--Take the privilege strings and convert them into privileges we can use to generate the task
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local read1_privs = terralib.newlist()
local write1_privs = terralib.newlist()
local update_neighbours = false
for _, v in pairs(read1) do
  read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
end
--Add the extra privilege required by the search algorithm
for _, v in pairs(write1) do
  write1_privs:insert( regentlib.privilege(regentlib.writes, parts1, string_to_field_path.get_field_path(v)))
  if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" then
    update_neighbours = true
  end
end
local task pairwise_task([parts1], config : region(ispace(int1d), config_type)) where
   [read1_privs], [write1_privs], reads(config) do
   for part1 in [parts1].ispace do
         [kernel_name(rexpr [parts1][part1] end, rexpr config[0] end)]
   end
end
return pairwise_task, update_neighbours
end

function run_per_particle_task( kernel_name )
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local cell1 = regentlib.newsymbol("cell1")
local cell_space = regentlib.newsymbol("cell_space")
local per_part = rquote

    --Do nothing
end

return run_per_particle_task_code
end


import "regent"

require("defaults")


--Generate the classic MD-style symmetric update pairwise task.
--This function assumes the cutoff is the same for both particles
function generate_symmetric_pairwise_task( kernel_name )

local task pairwise_task(parts1 : region(ispace(int1d),part), parts2 : region(ispace(int1d),part), space : region(ispace(int1d), space_config)) where 
   reads(parts1, parts2, space), writes(parts1, parts2) do
   for part1 in parts1 do
     for part2 in parts2 do
         [kernel_name(rexpr parts1[part1] end, rexpr parts2[part2] end) ]
     end
   end
end
return pairwise_task
end

--Functionality added but unclear on use-case right now
function generate_asymmetric_pairwise_task( kernel_name )
--Asymmetric kernel can only write to part1

local task pairwise_task(parts1 : region(ispace(int1d),part), parts2 : region(ispace(int1d),part), space : region(ispace(int1d), space_config)) where reads(parts1, parts2), writes(parts1) do
   for part1 in parts1 do
     for part2 in parts2 do
       [kernel_name( rexpr parts1[part1]end, rexpr parts2[part2] end)]
     end
   end
end
return pairwise_task
end

--TODO: Other pairwise tasks may be needed for other task-types (e.g. astro SPH)
--NYI

function run_symmetric_pairwise_task( task_name )
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local cell1 = regentlib.newsymbol("cell1")
local cell2 = regentlib.newsymbol("cell2")
local cell_space = regentlib.newsymbol("cell_space")
local task run_symmetric_pairwise_task_code( parts : region(ispace(int1d), part) )

    --Do nothing   
end

return run_symmetric_pairwise_task_code
end




-----------------------------------------
--End of Pair Tasks ---------------------
-----------------------------------------

--Generate a task to be executed on every particle in the system
function generate_per_part_task( kernel_name )

local task pairwise_task(parts1 : region(ispace(int1d),part), space : region(ispace(int1d), space_config)) where
   reads(parts1, space), writes(parts1) do
   for part1 in parts1 do
         [kernel_name(rexpr parts1[part1] end)]
   end
end
return pairwise_task
end

function run_per_particle_task( task_name )
local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local cell1 = regentlib.newsymbol("cell1")
local cell_space = regentlib.newsymbol("cell_space")
local task run_per_particle_task_code( parts : region(ispace(int1d), part) )

    --Do nothing
end

return run_per_particle_task_code
end


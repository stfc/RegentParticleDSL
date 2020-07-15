import "regent"

local format = require("std/format")
require("example_types/example_part")

--For best performance all tasks should use __demand(__inline) if possible.
--__demand(__inline)
--task example_pairwise_interaction(part1 : part, part2 : part)
--do

--We can access variables from the core_part or neighbour_part here
--  if(part1.core_part_space.pos_x > -1000.0) then
--    part1.extra_variable_3 = 1
--    part2.extra_variable_3 = 2
--  end
--
--end

--function example_pairwise_interaction( )
--
--local pairwise_interaction_function = terralib.newlist()
--
--pairwise_interaction_function:insert(rquote
--  if(part1.core_part_space.pos_x > -1000.0) then
--     part1.extra_variable_3 = 1
--     part2.extra_variable_3 = 2
--  end
--
--  end)
--
--  return pairwise_interaction_function
--end


function kernel_one(part1, part2)
local kernel = rquote
[part1].extra_variable_1 = 1.0
end
return kernel
end

function kernel_two(part1, part2)

local kernel = rquote
    [part2].extra_variable_1 = 2 
end 
return kernel
end

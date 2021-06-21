-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

require("defaults")


--Template for a symmetric pairwise kernel function
function pairwise_kernel(part1, part2, r2, config)
local kernel = rquote
    --Kernel code here, e.g.
    part1.extra_variable_1 = 1.0
    part2.extra_variable_1 = 1.0
end
return kernel
end

--Template for an asymmetric pairwise kernel function
function asym_pairwise_kernel(part1, part2, r2, config)
local kernel = rquote
    --Kernel code here, part2 is read-only,  e.g.
    part1.extra_variable_1 = part2.extra_variable_1 + 1.0
end
return kernel

--Template for a per-part kernel function
function per_part_kernel(part, config)

local kernel = rquote
    --Kernel code here, e.g.
    part.extra_variable_1 = 2
    part.core_part_space.vel_x *= 0.5 
end 
return kernel
end

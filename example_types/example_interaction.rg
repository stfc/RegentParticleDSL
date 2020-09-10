-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

local format = require("std/format")
require("example_types/example_part")


function kernel_one(part1, part2, r2)
local kernel = rquote
    part1.extra_variable_1 = 1.0
end
return kernel
end

function kernel_two(part1, part2, r2)

local kernel = rquote
    part2.extra_variable_1 = 2 
end 
return kernel
end

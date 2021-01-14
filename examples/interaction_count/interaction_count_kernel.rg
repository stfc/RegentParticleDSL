-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

--Symmetric interaction count kernel
function symmetric_interaction_count_kernel(part1, part2, r2)
local kernel = rquote
  part1.interactions = part1.interactions + 1
  part2.interactions = part2.interactions + 1
end
return kernel
end

--Asymmetric interaction count kernel.
function asymmetric_interaction_count_kernel(part1, part2, r2)
local kernel = rquote
  part1.interactions = part1.interactions + 1
end
return kernel
end

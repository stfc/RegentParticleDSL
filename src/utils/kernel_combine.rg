-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

local kernel_combine = {}

--Experimental function to test combining a set of kernels into a single kernel instead. Used automatically by the code.
function kernel_combine.combine_kernels( kernel_list )

   --Code that expects a kernel assumes the format is a function that returns an rquote, so we create that
   --by combining all the kernels in the kernel_list inside a new function, and return that function
    local function combine_kernels_innerfunc(part1, part2, r2)
        local combined_kernels_rquote = rquote
            [kernel_list:map( function(kernel)
                     return kernel(part1, part2, r2)
            end)];
        end
        return combined_kernels_rquote
    end

    return combined_kernels_innerfunc 
end

return kernel_combine

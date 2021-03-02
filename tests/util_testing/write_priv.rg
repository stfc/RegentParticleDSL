import "regent"

local compute_privileges = require("src/utils/compute_privilege")

local function symmetric_interaction_count_kernel(part1, part2, r2)
  local kernel = rquote
    part1.interaction = 2.5
    part2.interaction = 2.5
  end
return kernel
end

local function asymmetric_interaction_count_kernel(part1, part2, r2)
  local kernel = rquote
    part1.interaction = 2.5
  end
return kernel
end

local function test_symmetric_write()
  local read1, read2, write1, write2, reduc1, reduc2 = compute_privileges.two_region_privileges(symmetric_interaction_count_kernel)
  local count = 0
  for k,v in pairs(write1) do
    count = count + 1
    assert(v == "interaction")
  end
  assert(count == 1, "Write 1 did not contain 1 element")
  for k,v in pairs(write2) do
    count = count + 1
    assert(v == "interaction")
  end
  assert(count == 2, "Write 2 did not contain 1 element")
--Skip reads for now as write implies read in the current implementation
  for k, v in pairs(reduc1) do
    count = count + 1
  end
  assert(count == 2, "Reduc 1 not empty as expected")
  for k, v in pairs(reduc2) do
    count = count + 1
  end
  assert(count == 2, "Reduc 2 not empty as expected")
end

local function test_asymmetric_write()
  local read1, read2, write1, write2, reduc1, reduc2 = compute_privileges.two_region_privileges(asymmetric_interaction_count_kernel)
  local count = 0
  for k,v in pairs(write1) do
    count = count + 1
    assert(v == "interaction")
  end
  assert(count == 1, "Write 1 did not contain 1 element")
  for k,v in pairs(write2) do
    count = count + 1
  end
  assert(count == 1, "Write 2 not empty as expected")
--Skip reads for now as write implies read in the current implementation
  for k,v in pairs(read2) do
    count = count + 1
  end
  assert(count == 1, "Read 2 not empty as expected")
  for k, v in pairs(reduc1) do
    count = count + 1
  end
  assert(count == 1, "Reduc 1 not empty as expected")
  for k, v in pairs(reduc2) do
    count = count + 1
  end
  assert(count == 1, "Reduc 2 not empty as expected")
end

test_symmetric_write()
test_asymmetric_write()
print("Write privilege test successful")

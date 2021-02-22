import "regent"

local compute_privileges = require("src/utils/compute_privilege")

local function symmetric_interaction_count_kernel(part1, part2, r2)
  local kernel = rquote
    part1.interactions += 1
    part2.interactions += 1
  end
return kernel
end

local function asymmetric_interaction_count_kernel(part1, part2, r2)
  local kernel = rquote
    part1.interactions += 1
  end
return kernel
end


local function test_symmetric_reduc()
  local read1, read2, write1, write2, reduc1, reduc2 = compute_privileges.two_region_privileges(symmetric_interaction_count_kernel)
  local count = 0
  for k,v in pairs(read1) do
    count = count + 1
  end
  assert(count == 0, "Read 1 not empty as expected")
  for k,v in pairs(read2) do
    count = count + 1
  end
  assert(count == 0, "Read 2 not empty as expected")
  for k,v in pairs(write1) do
    count = count + 1
  end
  assert(count == 0, "Write 1 not empty as expected")
  for k,v in pairs(write2) do
    count = count + 1
  end
  assert(count == 0, "Write 2 not empty as expected")
  for k, v in pairs(reduc1) do
    count = count + 1
    assert(v[1] == "interactions")
    assert(v[2] == "+")
  end
  assert(count == 1, "Did not get 1 element in reduc1 as expected")
  for k, v in pairs(reduc2) do
    count = count + 1
    assert(v[1] == "interactions")
    assert(v[2] == "+")
  end
  assert(count == 2, "Did not get 1 element in reduc2 as expected")
end



local function test_asymmetric_reduc()
  local read1, read2, write1, write2, reduc1, reduc2 = compute_privileges.two_region_privileges(asymmetric_interaction_count_kernel)
  local count = 0
  for k,v in pairs(read1) do
    count = count + 1
  end
  assert(count == 0, "Read 1 not empty as expected")
  for k,v in pairs(read2) do
    count = count + 1
  end
  assert(count == 0, "Read 2 not empty as expected")
  for k,v in pairs(write1) do
    count = count + 1
  end
  assert(count == 0, "Write 1 not empty as expected")
  for k,v in pairs(write2) do
    count = count + 1
  end
  assert(count == 0, "Write 2 not empty as expected")
  for k, v in pairs(reduc1) do
    count = count + 1
    assert(v[1] == "interactions")
    assert(v[2] == "+")
  end
  assert(count == 1, "Did not get 1 element in reduc1 as expected")
  for k, v in pairs(reduc2) do
    count = count + 1
  end
  assert(count == 1, "Reduc2 not empty as expected")

end

test_symmetric_reduc()
test_asymmetric_reduc()
print("Test successful")

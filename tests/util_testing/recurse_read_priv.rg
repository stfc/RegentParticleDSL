import "regent"

local compute_privileges = require("src/utils/compute_privilege")

local function symmetric_interaction_count_kernel(part1, part2, r2)
  local kernel = rquote
    var a = part1.A.B.C.D.E.F
    var b = part2.A.B.C.D.E.F
  end
return kernel
end

local function test_symmetric_read()
  local read1, read2, write1, write2, reduc1, reduc2 = compute_privileges.two_region_privileges(symmetric_interaction_count_kernel)
  local count = 0
  for k,v in pairs(read1) do
    count = count + 1
    print(v)
    assert(v == "interaction")
  end
  assert(count == 1, "Read 1 did not contain 1 element")
  for k,v in pairs(read2) do
    count = count + 1
    assert(v == "interaction")
  end
  assert(count == 2, "Read 2 did not contain 1 element")
  for k,v in pairs(write1) do
    count = count + 1
  end
  assert(count == 2, "Write 1 not empty as expected")
  for k,v in pairs(write2) do
    count = count + 1
  end
  assert(count == 2, "Write 2 not empty as expected")
  for k, v in pairs(reduc1) do
    count = count + 1
  end
  assert(count == 2, "Reduc 1 not empty as expected")
  for k, v in pairs(reduc2) do
    count = count + 1
  end
  assert(count == 2, "Reduc 2 not empty as expected")
end

test_symmetric_read()
print("Recurse read fields successful")

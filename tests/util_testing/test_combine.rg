import "regent"

local compute_privileges = require("src/utils/compute_privilege")
local kernel_combine = require("src/utils/kernel_combine")

local function kernel_one(part1, part2, r2)
  local kernel = rquote
    var a = part1.interaction
  end
return kernel
end

local function kernel_two(part1, part2, r2)
  local kernel = rquote
    var a = part2.interaction
  end
return kernel

end

local function test_kernel_one_read()
  local read1, read2, write1, write2, reduc1, reduc2 = compute_privileges.two_region_privileges(kernel_one)
  local count = 0
  for k,v in pairs(read1) do
    count = count + 1
    assert(v == "interaction")
  end
  assert(count == 1, "Read 1 did not contain 1 element")
  for k,v in pairs(read2) do
    count = count + 1
  end
  assert(count == 1, "Read 2 did not contain 0 elements")
  for k,v in pairs(write1) do
    count = count + 1
  end
  assert(count == 1, "Write 1 not empty as expected")
  for k,v in pairs(write2) do
    count = count + 1
  end
  assert(count == 1, "Write 2 not empty as expected")
  for k, v in pairs(reduc1) do
    count = count + 1
  end
  assert(count == 1, "Reduc 1 not empty as expected")
  for k, v in pairs(reduc2) do
    count = count + 1
  end
  assert(count == 1, "Reduc 2 not empty as expected")
end

local function test_kernel_two_read()
  local read1, read2, write1, write2, reduc1, reduc2 = compute_privileges.two_region_privileges(kernel_two)
  local count = 0
  for k,v in pairs(read2) do
    count = count + 1
    assert(v == "interaction")
  end
  assert(count == 1, "Read 2 did not contain 1 element")
  for k,v in pairs(read1) do
    count = count + 1
  end
  assert(count == 1, "Read 1 did not contain 0 elements")
  for k,v in pairs(write1) do
    count = count + 1
  end
  assert(count == 1, "Write 1 not empty as expected")
  for k,v in pairs(write2) do
    count = count + 1
  end
  assert(count == 1, "Write 2 not empty as expected")
  for k, v in pairs(reduc1) do
    count = count + 1
  end
  assert(count == 1, "Reduc 1 not empty as expected")
  for k, v in pairs(reduc2) do
    count = count + 1
  end
  assert(count == 1, "Reduc 2 not empty as expected")
end

local function test_kernel_combined(combined)
  local read1, read2, write1, write2, reduc1, reduc2 = compute_privileges.two_region_privileges(combined)
  local count = 0
  for k,v in pairs(read2) do
    count = count + 1
    assert(v == "interaction")
  end
  assert(count == 1, "Read 2 did not contain 1 element")
  for k,v in pairs(read1) do
    count = count + 1
    assert(v == "interaction")
  end
  assert(count == 2, "Read 2 did not contain 1 elements")
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

test_kernel_one_read()
test_kernel_two_read()
kernels = {}
table.insert(kernels, kernel_one)
table.insert(kernels, kernel_two)
--kernels[1] = test_kernel_one_read()
--kernels[2] = test_kernel_two_read()
local combined = kernel_combine.combine_kernels(kernels)
local a = regentlib.newsymbol("a")
local b = regentlib.newsymbol("b")
local c = regentlib.newsymbol("c")
z = combined(a, b, c)
print(z)
test_kernel_combined(combined)
print("Combine kernel success")


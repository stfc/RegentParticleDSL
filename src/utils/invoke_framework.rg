-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

--Other headers
local compute_privileges = require("src/utils/compute_privilege")
--By default, asking for PAIRWISE gives a SYMMETRIC_PAIRWISE operation
SYMMETRIC_PAIRWISE = 1
PAIRWISE = 1
ASYMMETRIC_PAIRWISE = 2
PER_PART = 3

BARRIER = 100
NO_BARRIER = 101

MULTI_KERNEL=1000
SINGLE_KERNEL=1001

--FIXME: Handle reductions
local function is_safe_to_combine(kernel, combined_kernels) 
--For now we ignore this
  local pre_read1 = terralib.newlist()
  local pre_write1 = terralib.newlist()
  local hash_r1 = {}
  local hash_w1 = {}
  --Compute the read/write requirements for the already combined kernels
  for _, kernel in pairs(combined_kernels) do
    local temp_r1, temp_r2, temp_w1, temp_w2 = compute_privileges.two_region_privileges(kernel)
    --Merge the read/writes for this kernel with previous ones, keeping uniqueness
    for _,v in pairs(temp_r1) do
      if( not hash_r1[v]) then
        pre_read1:insert(v)
        hash_r1[v] = true
      end
    end
    for _,v in pairs(temp_r2) do
      if( not hash_r1[v]) then
        pre_read1:insert(v)
        hash_r1[v] = true
      end
    end
    for _,v in pairs(temp_w1) do
      if( not hash_w1[v]) then
        pre_write1:insert(v)
        hash_w1[v] = true
      end
    end
    for _,v in pairs(temp_w2) do
      if( not hash_w1[v]) then
        pre_write1:insert(v)
        hash_w1[v] = true
      end
    end
  end





local read1, read2, write1, write2 = compute_privileges.two_region_privileges( kernel )
local safe_to_combine = true
for _, v in pairs(write1) do
  if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" or v == "core_part_space.cutoff" then
    safe_to_combine = false
  end
  --Handle WaW or WaR dependencies
  if (hash_w1[v]) or (hash_r1[v]) then
    safe_to_combine = false
  end
end
for _, v in pairs(write2) do
  if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" or v == "core_part_space.cutoff" then
    safe_to_combine = false
  end
  --Handle WaW or WaR dependencies
  if (hash_w1[v]) or (hash_r1[v]) then
    safe_to_combine = false
  end
end
for _,v in pairs(read1) do
  --Handle RaW dependency
  if (hash_w1[v]) then
    safe_to_combine = false
  end
end
for _,v in pairs(read2) do
  --Handle RaW dependency
  if (hash_w1[v]) then
    safe_to_combine = false
  end
end
--TODO: Ignore this for now
safe_to_combine = false
return safe_to_combine
end

local function invoke_multikernel(config, ...)
  local end_barrier = true
  --Loop over the arguments and check for synchronicity (or other things to be added).
for i= 1, select("#",...) do
    if select(i, ...) == NO_BARRIER then
      end_barrier = false
    end
  end
  local kernels = {}
  local last_type = -1
  local quote_list = terralib.newlist()

  --Loop through the inputs and find all the tables.
  for i=1, select("#", ...) do
    local v = select(i, ...)
    if type(v) == "table" then
      if v[1] == nil or v[2] == nil then
        print("The arguments to invoke were incorrect")
        os.exit(1)
      end
      local func = v[1]
      local type_iterate = v[2]
 --I think we can refactor this using some functions to make the code cleaner, and just check if last_type == type_iterate. AC
      if type_iterate == SYMMETRIC_PAIRWISE then
        local safe_to_combine, can_be_combined = is_safe_to_combine(func, kernels)
        if safe_to_combine and last_type == SYMMETRIC_PAIRWISE then
          table.insert(kernels, func)
        elseif safe_to_combine then
          if #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
            quote_list:insert( create_asymmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == PER_PART then
            quote_list:insert( run_per_particle_task_multikernel( config, neighbour_init.cell_partition, unpack(kernels) ) )
          end

          kernels = {}
          table.insert(kernels, func)
          last_type = SYMMETRIC_PAIRWISE
        else --GENERATE PREVIOUS AND CURRENT
          if #kernels > 0 and last_type == SYMMETRIC_PAIRWISE then
            quote_list:insert( create_symmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
            quote_list:insert( create_asymmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == PER_PART then
            quote_list:insert( run_per_particle_task_multikernel( config, neighbour_init.cell_partition, unpack(kernels) ) )
          end
          if can_be_combined then
            kernels = {}
            table.insert(kernels, func)
            last_type = SYMMETRIC_PAIRWISE
          else
            quote_list:insert( create_symmetric_pairwise_runner(func, config, neighbour_init.cell_partition))
            kernels = {}
            last_type = -1
          end
        end
      elseif type_iterate == ASYMMETRIC_PAIRWISE then
        local safe_to_combine, can_be_combined = is_safe_to_combine(func, kernels)
        if safe_to_combine and last_type == ASYMMETRIC_PAIRWISE then
          table.insert(kernels, func)
        elseif safe_to_combine then
          if #kernels > 0 and last_type == SYMMETRIC_PAIRWISE then
            quote_list:insert( create_symmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == PER_PART then
            quote_list:insert( run_per_particle_task_multikernel( config, neighbour_init.cell_partition, unpack(kernels) ) )
          end

          kernels = {}
          table.insert(kernels, func)
          last_type = ASYMMETRIC_PAIRWISE
        else --GENERATE PREVIOUS AND CURRENT
          if #kernels > 0 and last_type == SYMMETRIC_PAIRWISE then
            quote_list:insert( create_symmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
            quote_list:insert( create_asymmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == PER_PART then
            quote_list:insert( run_per_particle_task_multikernel( config, neighbour_init.cell_partition, unpack(kernels) ) )
          end

          if can_be_combined then
            kernels = {}
            table.insert(kernels, func)
            last_type = ASYMMETRIC_PAIRWISE
          else
            quote_list:insert( create_asymmetric_pairwise_runner(func, config, neighbour_init.cell_partition))
            kernels = {}
            last_type = -1
          end
        end
      elseif type_iterate == PER_PART then
        local safe_to_combine, can_be_combined = is_safe_to_combine(func, kernels)
        if safe_to_combine and last_type == PER_PART then
          table.insert(kernels, func)
        elseif safe_to_combine then
          if #kernels > 0 and last_type == SYMMETRIC_PAIRWISE then
            quote_list:insert( create_symmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
            quote_list:insert( create_asymmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          end

          kernels = {}
          table.insert(kernels, func)
          last_type = PER_PART
        else --GENERATE PREVIOUS AND CURRENT
          if #kernels > 0 and last_type == SYMMETRIC_PAIRWISE then
            quote_list:insert( create_symmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
            quote_list:insert( create_asymmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
          elseif #kernels > 0 and last_type == PER_PART then
            quote_list:insert( run_per_particle_task_multikernel( config, neighbour_init.cell_partition, unpack(kernels) ) )
          end
          if can_be_combined then
            kernels = {}
            table.insert(kernels, func)
            last_type = PER_PART
          else
            quote_list:insert(  run_per_particle_task( func, config, neighbour_init.cell_partition ) )
            kernels = {}
            last_type = -1
          end
        end
      else
        print("The kernel type passed to invoke was not recognized")
        os.exit(1)
      end
    end
  end
  --Generate the final quote
  if #kernels > 0 and last_type == SYMMETRIC_PAIRWISE then
    quote_list:insert( create_symmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
  elseif #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
    quote_list:insert( create_asymmetric_pairwise_runner_multikernel(config, neighbour_init.cell_partition, unpack(kernels) ) )
  elseif #kernels > 0 and last_type == PER_PART then
    quote_list:insert( run_per_particle_task_multikernel( config, neighbour_init.cell_partition, unpack(kernels) ) )
  end
  local barrier_quote = rquote
  end
  if end_barrier then
    barrier_quote = rquote
      c.legion_runtime_issue_execution_fence(__runtime(), __context())
    end
  end
  local invoke_quote = rquote
    [quote_list];
    [barrier_quote];
  end
  return invoke_quote
end

local function invoke_per_kernel(config, ...)
  local end_barrier = true
  --Loop over the arguments and check for synchronicity (or other things to be added).
  for i= 1, select("#",...) do
    if select(i, ...) == NO_BARRIER then
      end_barrier = false
    end
  end

  local quote_list = terralib.newlist()
  --Loop through the inputs and find all the tables.
  for i= 1, select("#",...) do
     local v = select(i, ...)
    if type(v) == "table" then
      if v[1] == nil or v[2] == nil then
        print("The arguments to invoke were incorrect")
        os.exit(1)
      end
      local func = v[1]
      local type_iterate = v[2]
      if type_iterate == SYMMETRIC_PAIRWISE then
        quote_list:insert( create_symmetric_pairwise_runner(func, config, neighbour_init.cell_partition) )
      elseif type_iterate == ASYMMETRIC_PAIRWISE then
        quote_list:insert( create_asymmetric_pairwise_runner(func, config, neighbour_init.cell_partition) )
      elseif type_iterate == PER_PART then
        quote_list:insert( run_per_particle_task(func, config,  neighbour_init.cell_partition) )
      else
        print("The kernel type passed to invoke was not recognized")
        os.exit(1)
      end
    end
  end
  local barrier_quote = rquote

  end
  if end_barrier then
    barrier_quote = rquote
      c.legion_runtime_issue_execution_fence(__runtime(), __context())
    end
  end
  local invoke_quote = rquote
    [quote_list];
    [barrier_quote];
  end
return invoke_quote
end

function invoke(config, ...)
  local multi_kernel = true
  --Loop over the arguments and check for synchronicity (or other things to be added).
  for i= 1, select("#",...) do
    if select(i, ...) == SINGLE_KERNEL then
      multi_kernel = false
    end
  end
  if multi_kernel then
    return invoke_multikernel(config, ...)
  else
    return invoke_per_kernel(config, ...)
  end
end  

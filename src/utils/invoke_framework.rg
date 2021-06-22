-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

--Other headers
local c = regentlib.c
local compute_privileges = require("src/utils/compute_privilege")
local kernel_combine = require("src/utils/kernel_combine")
--By default, asking for PAIRWISE gives a SYMMETRIC_PAIRWISE operation
SYMMETRIC_PAIRWISE = 1
PAIRWISE = 1
ASYMMETRIC_PAIRWISE = 2
PER_PART = 3

BARRIER = 100
NO_BARRIER = 101

MULTI_KERNEL=1000
SINGLE_KERNEL=1001

local function is_safe_to_combine_perpart(kernel, combined_kernels)

    local hash_re2 = {}
    local reducs = {}
    local safe_to_combine = true
    --We only care about what happens in the config, particle-local operations are always ok to combine.
    for _, kernel in pairs(combined_kernels) do
        local temp_r1, temp_r2, temp_r3, temp_w1, temp_w2, temp_re1, temp_re2, temp_re3 = nil
        temp_r1, temp_r2, temp_w1, temp_w2, temp_re1, temp_re2 = compute_privileges.two_region_privileges(kernel)
        for _,v in pairs(temp_re2) do
           if( not hash_re1[v[1]]) then
               hash_re1[v[1]] = true
               table.insert(reducs, v)
           else -- if we hash is true we need to check its the same operator
               local found = false
               for _,v1 in pairs(reducs) do
                   if v1[1] == v[1] and v1[2] == v[2] then
                        found = true
                   elseif v1[1] == v[1] then
                        --Can't combine as we have a multi reduction operator on a single field already
                        safe_to_combine = false
                   end
               end
               if (not found) then
                   hash_re1[v[1]] = true
                   table.insert(reducs, v)
               end
           end
        end 

    end

    local read1, read2, read3, write1, write2, reduc1, reduc2, reduc3 = nil
    read1, read2, write1, write2, reduc1, reduc2 = compute_privileges.two_region_privileges( kernel )

    --Check read2 for RaRe conflicts
    for _, v in pairs(read2) do
        if hash_re2[v] then
            safe_to_combine = false
        end
    end
    --ReaRe conflicts are ok provided they have the same operator
    for _, v in pairs(reduc2) do
        if hash_re2[v[1]] then
            for _,v1 in pairs(reducs) do
                if v1[1] == v[1] and v1[2] ~= v[2] then
                    safe_to_combine = false
                end
            end
        end
    end

    --Can always be combined in theory
    local can_be_combined = true
    return safe_to_combine, can_be_combined
end

local function is_safe_to_combine_pairwise(kernel, combined_kernels)

  local hash_r1 = {}
  local hash_r3 = {}
  local hash_w1 = {}
  local hash_re1 = {}
  local hash_re3 = {}
  local reducs = {}
  local safe_to_combine = true

  --Compute the read/write requirements for the already combined kernels
    for _, kernel in pairs(combined_kernels) do
        local temp_r1, temp_r2, temp_r3, temp_w1, temp_w2, temp_re1, temp_re2, temp_re3 = nil
        temp_r1, temp_r2, temp_r3, temp_w1, temp_w2, temp_re1, temp_re2, temp_re3 = compute_privileges.three_region_privileges(kernel)

        --Merge the read/writes for this kernel with previous ones, keeping uniqueness
        for _,v in pairs(temp_r1) do
            if( not hash_r1[v]) then
                hash_r1[v] = true
            end
        end
        for _,v in pairs(temp_r2) do
            if( not hash_r1[v]) then
                hash_r1[v] = true
            end
        end
        for _,v in pairs(temp_r3) do
            if( not hash_r3[v]) then
                hash_r3[v] = true
            end
        end
        for _,v in pairs(temp_w1) do
            if( not hash_w1[v]) then
                hash_w1[v] = true
            end
        end
        for _,v in pairs(temp_w2) do
            if( not hash_w1[v]) then
                hash_w1[v] = true
            end
        end
        for _,v in pairs(temp_re1) do
            if( not hash_re1[v[1]]) then
                hash_re1[v[1]] = true
                table.insert(reducs, v)
            else -- if we hash is true we need to check its the same operator
                local found = false
                for _,v1 in pairs(reducs) do
                    if v1[1] == v[1] and v1[2] == v[2] then
                      found = true
                    elseif v1[1] == v[1] then
                      --Can't combine as we have a multi reduction operator on a single field already
                      safe_to_combine = false
                    end
                end
                if (not found) then
                    hash_re1[v[1]] = true
                    table.insert(reducs, v)
                end
            end
        end
        for _,v in pairs(temp_re2) do
           if( not hash_re1[v[1]]) then
               hash_re1[v[1]] = true
               table.insert(reducs, v)
           else -- if we hash is true we need to check its the same operator
               local found = false
               for _,v1 in pairs(reducs) do
                   if v1[1] == v[1] and v1[2] == v[2] then
                     found = true
                   elseif v1[1] == v[1] then
                     --Can't combine as we have a multi reduction operator on a single field already
                     safe_to_combine = false
                   end
               end
               if (not found) then
                   hash_re1[v[1]] = true
                   table.insert(reducs, v)
               end
           end
        end
        for _,v in pairs(temp_re3) do
           if( not hash_re3[v[1]]) then
               hash_re3[v[1]] = true
               table.insert(reducs, v)
           else -- if we hash is true we need to check its the same operator
               local found = false
               for _,v1 in pairs(reducs) do
                   if v1[1] == v[1] and v1[2] == v[2] then
                     found = true
                   elseif v1[1] == v[1] then
                     --Can't combine as we have a multi reduction operator on a single field already
                     safe_to_combine = false
                   end
               end
               if (not found) then
                   hash_re3[v[1]] = true
                   table.insert(reducs, v)
               end
           end
        end
    end


    local read1, read2, read3, write1, write2, reduc1, reduc2, reduc3 = nil
    read1, read2, read3, write1, write2, reduc1, reduc2, reduc3 =  compute_privileges.three_region_privileges(kernel)

    local can_be_combined = true
    for _, v in pairs(write1) do
      if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" or v == "core_part_space.cutoff" then
        safe_to_combine = false
        can_be_combined = false
      end
      --Handle WaW or WaR dependencies
      if (hash_w1[v]) or (hash_r1[v]) or (hash_re1[v]) then
        safe_to_combine = false
      end
    end
    for _, v in pairs(write2) do
      if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" or v == "core_part_space.cutoff" then
        safe_to_combine = false
        can_be_combined = false
      end
      --Handle WaW or WaR dependencies
      if (hash_w1[v]) or (hash_r1[v]) or (hash_re1[v]) then
        safe_to_combine = false
      end
    end
    for _,v in pairs(read1) do
      --Handle RaW dependency
      if (hash_w1[v]) or (hash_re1[v]) then
        safe_to_combine = false
      end
    end
    for _,v in pairs(read2) do
      --Handle RaW dependency
      if (hash_w1[v]) or (hash_re1[v]) then
        safe_to_combine = false
      end
    end
    for _, v in pairs(read3) do
      --Handle RaW dependency
      if (hash_re3[v]) then
        safe_to_combine = false
      end
    end

    --Check the reduction operations are ok
    for _, v in pairs(reduc1) do
      if v[1] =="core_part_space.pos_x" or v[1] == "core_part_space.pos_y" or v[1] == "core_part_space.pos_z" or v[1] == "core_part_space.cutoff" then
        safe_to_combine = false
        can_be_combined = false
      end
      --Check if ReaW or ReaR
      if (hash_w1[v[1]]) or (hash_r1[v[1]]) then
        safe_to_combine = false
      end
      --Check reduction compatability
      if (hash_re1[v[1]]) then
          --Check if the operator is the same
          for _,v1 in pairs(reducs) do
              if v1[1] == v[1] and v1[2] ~= v[2] then
                safe_to_combine = false
              end
              --We also need to add this into check (since another reduction on the same thing would hurt parallelism if we combine)
              if v1[1] ~= v[1] then
                hash_re1[v[1]] = true
                table.insert(reducs, v)
              end
          end
      end
    end
    for _, v in pairs(reduc2) do
      if v[1] =="core_part_space.pos_x" or v[1] == "core_part_space.pos_y" or v[1] == "core_part_space.pos_z" or v[1] == "core_part_space.cutoff" then
        safe_to_combine = false
        can_be_combined = false
      end
      --Check if ReaW or ReaR
      if (hash_w1[v[1]]) or (hash_r1[v[1]]) then
        safe_to_combine = false
      end
      --Check reduction compatability
      if (hash_re1[v[1]]) then
          --Check if the operator is the same
          for _,v1 in pairs(reducs) do
              if v1[1] == v[1] and v1[2] ~= v[2] then
                safe_to_combine = false
              end
              --We also need to add this into check (since another reduction on the same thing would hurt parallelism if we combine)
              if v1[1] ~= v[1] then
                hash_re1[v[1]] = true
                table.insert(reducs, v)
              end
          end
      end
    end
    --Only per-part tasks can be combined if they affect the config
    for _,v in pairs(reduc3) do
        safe_to_combine = false
        can_be_combined = false
    end

    return safe_to_combine, can_be_combined
end

--FIXME: Handle reductions
local function is_safe_to_combine(kernel, combined_kernels, type_iterate) 
    local safe_to_combine  = false
    local can_be_combined = false
    if type_iterate == PER_PART then
        safe_to_combine, can_be_combined = is_safe_to_combine_perpart(kernel, combined_kernels)
    elseif type_iterate == SYMMETRIC_PAIRWISE or type_iterate == ASYMMETRIC_PAIRWISE then
        safe_to_combine, can_be_combined = is_safe_to_combine_pairwise(kernel, combined_kernels)
    end

return safe_to_combine, can_be_combined
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
        local safe_to_combine, can_be_combined = is_safe_to_combine(func, kernels, type_iterate)
        if safe_to_combine and last_type == SYMMETRIC_PAIRWISE then
          table.insert(kernels, func)
        elseif safe_to_combine then
          if #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
            local combined_kernel = kernel_combine.combine_kernels(kernels)
            quote_list:insert( create_assymetric_pairwise_runner(combined_kernel, config, neighbour_init.cell_partition) )
          elseif #kernels > 0 and last_type == PER_PART then
            local combined_kernel = kernel_combine.combine_kernels(kernels)
            quote_list:insert( run_per_particle_task( combined_kernel, config, neighbour_init.cell_partition ) )
          end

          kernels = {}
          table.insert(kernels, func)
          last_type = SYMMETRIC_PAIRWISE
        else --GENERATE PREVIOUS AND CURRENT
          if #kernels > 0 and last_type == SYMMETRIC_PAIRWISE then
            local combined_kernel = kernel_combine.combine_kernels(kernels)
            quote_list:insert( create_symmetric_pairwise_runner( combined_kernel, config, neighbour_init.cell_partition ) )
          elseif #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
            local combined_kernel = kernel_combine.combine_kernels(kernels)
            quote_list:insert( create_asymmetric_pairwise_runner( combined_kernel, config, neighbour_init.cell_partition) )
          elseif #kernels > 0 and last_type == PER_PART then
            local combined_kernel = kernel_combine.combine_kernels(kernels)
            quote_list:insert( run_per_particle_task( combined_kernel, config, neighbour_init.cell_partition) )
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
        local safe_to_combine, can_be_combined = is_safe_to_combine(func, kernels, type_iterate)
        if safe_to_combine and last_type == ASYMMETRIC_PAIRWISE then
          table.insert(kernels, func)
        elseif safe_to_combine then
          if #kernels > 0 and last_type == SYMMETRIC_PAIRWISE then
            local combined_kernel = kernel_combine.combine_kernels(kernels)
            quote_list:insert( create_symmetric_pairwise_runner( combined_kernel, config, neighbour_init.cell_partition ) ) 
          elseif #kernels > 0 and last_type == PER_PART then
            local combined_kernel = kernel_combine.combine_kernels(kernels)
            quote_list:insert( run_per_particle_task( combined_kernel, config, neighbour_init.cell_partition ) ) 
          end

          kernels = {}
          table.insert(kernels, func)
          last_type = ASYMMETRIC_PAIRWISE
        else --GENERATE PREVIOUS AND CURRENT
          if #kernels > 0 and last_type == SYMMETRIC_PAIRWISE then
            local combined_kernel = kernel_combine.combine_kernels(kernels)
            quote_list:insert( create_symmetric_pairwise_runner( combined_kernel, config, neighbour_init.cell_partition ) )
          elseif #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
            local combined_kernel = kernel_combine.combine_kernels(kernels)
            quote_list:insert( create_asymmetric_pairwise_runner(combined_kernel, config, neighbour_init.cell_partition ) )
          elseif #kernels > 0 and last_type == PER_PART then
            local combined_kernel = kernel_combine.combine_kernels(kernels)
            quote_list:insert( run_per_particle_task( combined_kernel, config, neighbour_init.cell_partition ) )
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
        local safe_to_combine, can_be_combined = is_safe_to_combine(func, kernels, type_iterate)
        if safe_to_combine and last_type == PER_PART then
          table.insert(kernels, func)
        elseif safe_to_combine then
          if #kernels > 0 and last_type == SYMMETRIC_PAIRWISE then
            local combined_kernel = kernel_combine.combine_kernels(kernels)
            quote_list:insert( create_symmetric_pairwise_runner(combined_kernel, config, neighbour_init.cell_partition) )
          elseif #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
            local combined_kernel = kernel_combine.combine_kernels(kernels)
            quote_list:insert( create_asymmetric_pairwise_runner(combined_kernel, config, neighbour_init.cell_partition) )
          end

          kernels = {}
          table.insert(kernels, func)
          last_type = PER_PART
        else --GENERATE PREVIOUS AND CURRENT
          if #kernels > 0 and last_type == SYMMETRIC_PAIRWISE then
            local combined_kernel = kernel_combine.combine_kernels(kernels)
            quote_list:insert( create_symmetric_pairwise_runner(combined_kernel, config, neighbour_init.cell_partition) )
          elseif #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
            local combined_kernel = kernel_combine.combine_kernels(kernels)
            quote_list:insert( create_asymmetric_pairwise_runner(combined_kernel, config, neighbour_init.cell_partition ) )
          elseif #kernels > 0 and last_type == PER_PART then
            local combined_kernel = kernel_combine.combine_kernels(kernels)
            quote_list:insert( run_per_particle_task( combined_kernel, config, neighbour_init.cell_partition) )
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
    local combined_kernel = kernel_combine.combine_kernels(kernels)
    quote_list:insert( create_symmetric_pairwise_runner(combined_kernel, config, neighbour_init.cell_partition ) )
  elseif #kernels > 0 and last_type == ASYMMETRIC_PAIRWISE then
    local combined_kernel = kernel_combine.combine_kernels(kernels)
    quote_list:insert( create_asymmetric_pairwise_runner(combined_kernel, config, neighbour_init.cell_partition ) )
  elseif #kernels > 0 and last_type == PER_PART then
    local combined_kernel = kernel_combine.combine_kernels(kernels)
    quote_list:insert( run_per_particle_task( combined_kernel, config, neighbour_init.cell_partition ) )
  end
  local barrier_quote = rquote
  end
  if DSL_settings.TIMING or end_barrier then
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
  if DSL_settings.TIMING or end_barrier then
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

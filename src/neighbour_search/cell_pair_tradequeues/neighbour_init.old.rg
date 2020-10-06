import "regent"

require("defaults")
require("src/neighbour_search/cell_pair_tradequeues/cell")

local string_to_field_path = require("src/utils/string_to_fieldpath")
local recursive_fields = require("src/utils/recursive_fields")
local neighbour_init = {}

neighbour_init.padded_particle_array = regentlib.newsymbol("padded_particle_array")
neighbour_init.TradeQueueRegion = regentlib.newsymbol("TradeQueueRegion")
neighbour_init.TradeQueues = regentlib.newsymbol("TradeQueues")
neighbour_init.cell_partition = regentlib.newsymbol("padded_cell_partition")

local DEBUG = true

local directions = {} --Lookup table for the directions
directions[0] =  {-1, -1, -1}
directions[1] =  {-1, -1,  0}
directions[2] =  {-1, -1,  1}
directions[3] =  {-1,  0, -1}
directions[4] =  {-1,  0,  0}
directions[5] =  {-1,  0,  1}
directions[6] =  {-1,  1, -1}
directions[7] =  {-1,  1,  0}
directions[8] =  {-1,  1,  1}
directions[9] =  { 0, -1, -1}
directions[10] = { 0, -1,  0}
directions[11] = { 0, -1,  1}
directions[12] = { 0,  0, -1} --We skip 0,0,0 as cell with itself is not relevant for tradequeues
directions[13] = { 0,  0,  1}
directions[14] = { 0,  1, -1}
directions[15] = { 0,  1,  0}
directions[16] = { 0,  1,  1}
directions[17] = { 1, -1, -1}
directions[18] = { 1, -1,  0}
directions[19] = { 1, -1,  1}
directions[20] = { 1,  0, -1}
directions[21] = { 1,  0,  0}
directions[22] = { 1,  0,  1}
directions[23] = { 1,  1, -1}
directions[24] = { 1,  1,  0}
directions[25] = { 1,  1,  1}

--Must be a better way to implement this...
__demand(__inline)
task lookup_shift(direction : int) : int3d
  var rval = int3d({0,0,0})
  if direction == 0 then
    rval = int3d({-1, -1, -1})
  elseif direction == 1 then
    rval = int3d({-1, -1,  0})
  elseif direction == 2 then
    rval = int3d({-1, -1,  1})
  elseif direction == 3 then
    rval = int3d({-1,  0, -1})
  elseif direction == 4 then
    rval = int3d({-1,  0,  0})
  elseif direction == 5 then
    rval = int3d({-1,  0,  1})
  elseif direction == 6 then
    rval = int3d({-1,  1, -1})
  elseif direction == 7 then
    rval = int3d({-1,  1,  0})
  elseif direction == 8 then
    rval = int3d({-1,  1,  1})
  elseif direction == 9 then
    rval = int3d({ 0, -1, -1})
  elseif direction == 10 then
    rval = int3d({ 0, -1,  0})
  elseif direction == 11 then
    rval = int3d({ 0, -1,  1})
  elseif direction == 12 then
    rval = int3d({ 0,  0, -1})
  elseif direction == 13 then
    rval = int3d({ 0,  0,  1})
  elseif direction == 14 then
    rval = int3d({ 0,  1, -1})
  elseif direction == 15 then
    rval = int3d({ 0,  1,  0})
  elseif direction == 16 then
    rval = int3d({ 0,  1,  1})
  elseif direction == 17 then
    rval = int3d({ 1, -1, -1})
  elseif direction == 18 then
    rval = int3d({ 1, -1,  0})
  elseif direction == 19 then
    rval = int3d({ 1, -1,  1})
  elseif direction == 20 then
    rval = int3d({ 1,  0, -1})
  elseif direction == 21 then
    rval = int3d({ 1,  0,  0})
  elseif direction == 22 then
    rval = int3d({ 1,  0,  1})
  elseif direction == 23 then
    rval = int3d({ 1,  1, -1})
  elseif direction == 24 then
    rval = int3d({ 1,  1,  0})
  elseif direction == 25 then
    rval = int3d({ 1,  1,  1})
  end
  return rval
end


local function construct_part_structure()
  local part_structure = terralib.newlist()
  local field_strings = {}
  local type_table = {}
  for k, v in pairs(part.fields) do
    recursive_fields.recurse_field(v, field_strings, type_table)
  end
  for k, _ in pairs(field_strings) do
    part_structure:insert({field = string_to_field_path.get_field_path(field_strings[k])})
  end
  return part_structure
end
local part_structure = construct_part_structure()

--FIXME: Optimisation of this value might be good, for now we're allowing cells to grow by 50 particles before needing a repartition.
--This realistically is simulation dependent
neighbour_init.padding_per_cell = 50
neighbour_init.tradequeue_size = neighbour_init.padding_per_cell / 5

--This function updates the cells to reflect any motion that occurs in the system. We repartition only as required, but never change
--the number of cells at this point due to causing issues with various assumptions in the system at the moment.
function neighbour_init.update_cells(variables)

local task squash_cell(cell : region(ispace(int1d), part) ) where reads(cell), writes(cell) do

  var lo = cell.ispace.bounds.lo
  var hi = cell.ispace.bounds.hi
  var last_valid = hi

  for x= int32(hi), int32(lo)-1, -1 do
    if(not cell[int1d(x)].neighbour_part_space._valid) then
      --If we thought this was valid but isn't, then move the valid counter backwards
      if int1d(x) == last_valid then
        last_valid = last_valid - int1d(1)
      else
      --Otherwise swap it with the last valid entry
        [part_structure:map(function(element)
          return rquote
            cell[x].[element.field] = cell[last_valid].[element.field]
          end
        end)];
        --Set the last valid we copied to invalid
        cell[last_valid].neighbour_part_space._valid = false
        var changed = false
        -- Work backwards to find the new last valid, we know that x is always valid
          for z = int32(last_valid-1), x-1, -1 do
            --Once we find then we don't continue to search - no break statement in Lua - avoiding using goto
            if(not changed) then
              if cell[int1d(z)].neighbour_part_space._valid then
                last_valid = int1d(z)
                changed = true
              end
            end
          end
      end
    end
  end
end

--__demand(__leaf)
local task compute_new_dests(particles : region(ispace(int1d), part), config : region(ispace(int1d), config_type)) where
  reads(particles, config), writes(particles.neighbour_part_space.cell_id, particles.core_part_space) do
  for particle in particles do
    --Ignore non-valid particles
    if (particles[particle].neighbour_part_space._valid) then
      if particles[particle].core_part_space.id == int1d(8) then
        format.println("Part 8 in {} {} {}", particles[particle].neighbour_part_space.cell_id.x, 
                                             particles[particle].neighbour_part_space.cell_id.y, 
                                             particles[particle].neighbour_part_space.cell_id.z)
      end
          --We need to handle periodicity
      if (particles[particle].core_part_space.pos_x >= config[0].space.dim_x) then 
        particles[particle].core_part_space.pos_x = particles[particle].core_part_space.pos_x - config[0].space.dim_x
      end
      if (particles[particle].core_part_space.pos_y >= config[0].space.dim_y) then 
        particles[particle].core_part_space.pos_y = particles[particle].core_part_space.pos_y - config[0].space.dim_y
      end
      if (particles[particle].core_part_space.pos_z >= config[0].space.dim_z) then 
        particles[particle].core_part_space.pos_z = particles[particle].core_part_space.pos_z - config[0].space.dim_z
      end
      if (particles[particle].core_part_space.pos_x < 0.0) then 
        particles[particle].core_part_space.pos_x = particles[particle].core_part_space.pos_x + config[0].space.dim_x
      end
      if (particles[particle].core_part_space.pos_y < 0.0) then 
        particles[particle].core_part_space.pos_y = particles[particle].core_part_space.pos_y + config[0].space.dim_y
      end
      if (particles[particle].core_part_space.pos_z < 0.0) then 
        particles[particle].core_part_space.pos_z = particles[particle].core_part_space.pos_z + config[0].space.dim_z
      end
      var x_cell : int1d = int1d( (particles[particle].core_part_space.pos_x / config[0].neighbour_config.cell_dim_x))
      var y_cell : int1d = int1d( (particles[particle].core_part_space.pos_y / config[0].neighbour_config.cell_dim_y))
      var z_cell : int1d = int1d( (particles[particle].core_part_space.pos_z / config[0].neighbour_config.cell_dim_z))
      var cell_loc : int3d = int3d( {x_cell, y_cell, z_cell} )
      particles[particle].neighbour_part_space.cell_id = cell_loc
      if particles[particle].core_part_space.id == int1d(8) then
        format.println("Part 8 in {} {} {}", x_cell, y_cell, z_cell)
      end
    end
  end
end

--__demand(__leaf)
local task tradequeue_push(parts : region(ispace(int1d), part), tradequeue : region(ispace(int1d), part), cell_id : int3d, neighbour : int3d) : bool where
    reads(parts, tradequeue), writes(tradequeue, parts) do

    --Keep track of whether the tradequeues succeed
    var valid = true
    --Find end of array
    var hi = int32(parts.ispace.bounds.hi)
    --Find start of array
    var lo = int32(parts.ispace.bounds.lo)
    var last_valid = parts.ispace.bounds.hi 
    var tradequeue_added = 0
--    format.println(" {} {} {}", hi, lo, last_valid)
--    if cell_id == int3d({2,2,0}) then
--      format.println("cell_id {} {} {}, neighbour {} {} {}", cell_id.x, cell_id.y, cell_id.z, neighbour.x, neighbour.y, neighbour.z)
--    end
      for i in parts.ispace do
        if parts[i].neighbour_part_space._valid and parts[i].core_part_space.id == int1d(8) then
          format.println("PUSH FOUND PART 8 IN CELL {} {} {} ({} {})",cell_id.x, cell_id.y, cell_id.z, lo, hi)
        elseif not(parts[i].neighbour_part_space._valid) and parts[i].core_part_space.id == int1d(8) then
          format.println("PUSH INVALID PART 8 IN CELL {} {} {}",cell_id.x, cell_id.y, cell_id.z)
        end
      end
    --Inside terra so C style loop, x=hi, x > lo, x--
    for part=hi, lo-1, -1 do
--      format.println("{}", part)
      if(parts[int1d(part)].neighbour_part_space._valid) then
        --If we're moving this particle to the trade queue
--        format.println("part {} {} {} neighbour {} {} {}", parts[int1d(part)].neighbour_part_space.cell_id.x,
--                                                           parts[int1d(part)].neighbour_part_space.cell_id.y,
--                                                           parts[int1d(part)].neighbour_part_space.cell_id.z,
--                                                           neighbour.x,
--                                                           neighbour.y,
--                                                           neighbour.z)
        if parts[int1d(part)].neighbour_part_space.cell_id == neighbour then
          if parts[int1d(part)].core_part_space.id == int1d(8)  then
             format.println("moving particle {}, validity {} to cell {} {} {} from {} {} {} at index: {} of ({} {})", parts[int1d(part)].core_part_space.id,
                                                                             [int32](parts[int1d(part)].neighbour_part_space._valid),
                                                                             neighbour.x, neighbour.y, neighbour.z, cell_id.x, cell_id.y, cell_id.z, part, lo, hi)
           end
          if(tradequeue_added < neighbour_init.tradequeue_size) then
            --Copy data to the trade queue
            [part_structure:map(function(element)
              return rquote
                tradequeue[int1d(tradequeue_added)].[element.field] = parts[int1d(part)].[element.field]
              end
            end)];
            tradequeue_added = tradequeue_added + 1
            --If we are the last valid particle we only decrement the last_valid instead of also
            --copying it.
            if(last_valid ~= int1d(part)) then
              --Copy last_valid to this index.
              [part_structure:map(function(element)
                return rquote
                parts[int1d(part)].[element.field] = parts[last_valid].[element.field]
                end
              end)];
              parts[last_valid].neighbour_part_space._valid = false
              last_valid = last_valid - int1d(1)
            else
             --Mark as invalid and reduce counter
             parts[int1d(part)].neighbour_part_space._valid = false
             last_valid = int1d(part-1)
            end
          else
            valid = false
          end
        end
      else
        last_valid = int1d(part-1)  
      end
    end
    return valid
end

--If we're running with debugging then we add an assertion to check for contiguity
local contig_assertion = function()
return rquote end
end
if DEBUG ~= nil then
contig_assertion = function(first_invalid, part, parts, cell)
local rval = rquote
  if(first_invalid ~= part) then
    format.println("first_invalid {}, part {}", first_invalid, part)
    for x = 0, first_invalid+1, 1 do
      format.println("Validity of index {} is {}, has id {}", x, [int32](parts[int1d(x)].neighbour_part_space._valid), 
                                                                      parts[int1d(x)].core_part_space.id) 
    end
    format.println("Part at position {} has ID {}", part, parts[int1d(part)].core_part_space.id)
    format.println("Particles in {} {} {} are not contiguous", cell.x,cell.y, cell.z)
    regentlib.assert(first_invalid == part, "The particles in this cell are not contiguous which is unexpected.")
  end
  end
  return rval
end
end 

--__demand(__leaf)
local task tradequeue_pull(parts : region(ispace(int1d), part), tradequeue : region(ispace(int1d), part), cell : int3d) : bool where
    reads(parts, tradequeue), writes(tradequeue, parts) do

    --Keep track of whether the tradequeues succeed
    var valid = true
    --Find end of array
    var hi = int32(parts.ispace.bounds.hi)
    --Find start of array
    var lo = int32(parts.ispace.bounds.lo)
    var first_invalid = lo
    --Find the first index containing no valid particle in the cell
    __forbid(__vectorize)
    for part = lo, hi+1 do
      if parts[int1d(part)].neighbour_part_space._valid then
        [contig_assertion(first_invalid, part, parts, cell)];
        first_invalid = part + 1
      end 
    end
    for part in tradequeue.ispace do
      if valid and tradequeue[part].neighbour_part_space._valid and tradequeue[part].neighbour_part_space.cell_id == cell then
        if first_invalid <= hi then
          --Copy the particle data into the first invalid index
            [part_structure:map(function(element)
              return rquote
                parts[first_invalid].[element.field] = tradequeue[part].[element.field]
              end
            end)];
          format.println("Copying {} into cell {} {} {} at index {}", parts[first_invalid].core_part_space.id, 
                                                                      cell.x, cell.y, cell.z,
                                                                      first_invalid)
          first_invalid = first_invalid + int1d(1)
          tradequeue[part].neighbour_part_space._valid = false
        else
          valid = false
        end
      end
    end

    return valid
end

local assert_correct_cells = function()
  return rquote
  end
end
if DEBUG ~= nil then
  assert_correct_cells = function()
    local rval = rquote
      --format.println("Checking for correctness")
      for cell in neighbour_init.cell_partition.colors do
        for part in neighbour_init.cell_partition[cell] do
          if neighbour_init.cell_partition[cell][part].neighbour_part_space._valid then
            if neighbour_init.cell_partition[cell][part].neighbour_part_space.cell_id ~= int3d(cell) then
               format.println("Particle {} meant to be in cell {} {} {} but in cell {} {} {}", 
                               neighbour_init.cell_partition[cell][part].core_part_space.id,
                               neighbour_init.cell_partition[cell][part].neighbour_part_space.cell_id.x,
                               neighbour_init.cell_partition[cell][part].neighbour_part_space.cell_id.y,
                               neighbour_init.cell_partition[cell][part].neighbour_part_space.cell_id.z,
                               cell.x, cell.y, cell.z)
               format.println("Position {} {} {}", neighbour_init.cell_partition[cell][part].core_part_space.pos_x,
                                                   neighbour_init.cell_partition[cell][part].core_part_space.pos_y,
                                                   neighbour_init.cell_partition[cell][part].core_part_space.pos_z)

            end
            regentlib.assert(neighbour_init.cell_partition[cell][part].neighbour_part_space.cell_id == int3d(cell), "particle found in wrong cell")
          end
        end
      end    
    end
    return rval
  end
end


local update_cells_quote = rquote

  var valid : bool = true
  for cell in [neighbour_init.cell_partition].colors do
    for i in [neighbour_init.cell_partition][cell].ispace do
      if [neighbour_init.cell_partition][cell][i].neighbour_part_space._valid and [neighbour_init.cell_partition][cell][i].core_part_space.id == int1d(8) then
        format.println("PART 8 IN UPDATECELL {} {} {}",cell.x, cell.y, cell.z)
      end
    end
  end
  __demand(__index_launch)
  for cell in neighbour_init.cell_partition.colors do
    compute_new_dests([neighbour_init.cell_partition][cell], [variables.config])
  end
  for direction = 0, 26 do
    if(valid) then
      var shift : int3d = lookup_shift(direction)
      --Push and pull tradequeues
      --Compute the neighbour in the relevant direction
      --First spawn all the push tasks
      for cell in [neighbour_init.cell_partition].colors do
        var new_x = cell.x + shift.x 
        var new_y = cell.y + shift.y
        var new_z = cell.z + shift.z
        if new_x < 0 then new_x = new_x + [variables.config][0].neighbour_config.x_cells end
        if new_y < 0 then new_y = new_y + [variables.config][0].neighbour_config.y_cells end
        if new_z < 0 then new_z = new_z + [variables.config][0].neighbour_config.z_cells end
        if new_x >= [variables.config][0].neighbour_config.x_cells then new_x = new_x - [variables.config][0].neighbour_config.x_cells end
        if new_y >= [variables.config][0].neighbour_config.y_cells then new_y = new_y - [variables.config][0].neighbour_config.y_cells end
        if new_z >= [variables.config][0].neighbour_config.z_cells then new_z = new_z - [variables.config][0].neighbour_config.z_cells end
        var neighbour : int3d = int3d({new_x, new_y, new_z})
--        if cell == int3d({2, 2, 0}) then
--          format.println("Cell {} {} {} neighbour {} {} {}", cell.x, cell.y, cell.z, neighbour.x, neighbour.y, neighbour.z)
--        end
        var tradequeueid = cell.x + cell.y * [variables.config][0].neighbour_config.x_cells + cell.z * [variables.config][0].neighbour_config.x_cells * [variables.config][0].neighbour_config.y_cells
        valid = tradequeue_push([neighbour_init.cell_partition][cell], neighbour_init.TradeQueues[tradequeueid], cell, neighbour) and valid
      end
      --Once all the push tasks are spawned, spawn the pull tasks
      for cell in [neighbour_init.cell_partition].colors do
        var new_x = cell.x + shift.x 
        var new_y = cell.y + shift.y
        var new_z = cell.z + shift.z
        if new_x < 0 then new_x = new_x + [variables.config][0].neighbour_config.x_cells end
        if new_y < 0 then new_y = new_y + [variables.config][0].neighbour_config.y_cells end
        if new_z < 0 then new_z = new_z + [variables.config][0].neighbour_config.z_cells end
        if new_x >= [variables.config][0].neighbour_config.x_cells then new_x = new_x - [variables.config][0].neighbour_config.x_cells end
        if new_y >= [variables.config][0].neighbour_config.y_cells then new_y = new_y - [variables.config][0].neighbour_config.y_cells end
        if new_z >= [variables.config][0].neighbour_config.z_cells then new_z = new_z - [variables.config][0].neighbour_config.z_cells end
        var neighbour : int3d = int3d({new_x, new_y, new_z})
        var tradequeueid = cell.x + cell.y * [variables.config][0].neighbour_config.x_cells + cell.z * [variables.config][0].neighbour_config.x_cells * [variables.config][0].neighbour_config.y_cells
        valid = valid and tradequeue_pull([neighbour_init.cell_partition][neighbour], neighbour_init.TradeQueues[tradequeueid], neighbour)
      end
    end
  end
  --Trade queues attempted
  if(valid) then
    [assert_correct_cells()];
  else
  --TODO: This may need serious optimisation for occasions where it is called regularly, however we hope this never happens "too often".
  --If this is a major issue for your code, please let us know on the issue:
  --https://github.com/stfc/RegentParticleDSL/issues/56
  --Make the particle array contiguous - we just use the squash_cell for now
    squash_cell([neighbour_init.padded_particle_array])
  --Find the end of the particles
  --If we have debug this will check for continuity
  --Find end of array
    var hi = int32([neighbour_init.padded_particle_array].ispace.bounds.hi)
    --Find start of array
    var lo = int32([neighbour_init.padded_particle_array].ispace.bounds.lo)
    var start_index = lo
    var temp :int3d = int3d({0,0,0})
    --Find the first index containing no valid particle in the cell
    __forbid(__vectorize)
    for part = lo, hi+1 do
      if [neighbour_init.padded_particle_array][int1d(part)].neighbour_part_space._valid then
        [contig_assertion(start_index, part, neighbour_init.padded_particle_array, temp)];
        start_index = part + 1
      end
    end
    --Now loop over the trade queues to find any lost particles
    for part in [neighbour_init.TradeQueueRegion].ispace do
      --If its valid, then we need to add it into the particle array
      if [neighbour_init.TradeQueueRegion][part].neighbour_part_space._valid then
           [part_structure:map(function(element)
              return rquote
                [neighbour_init.padded_particle_array][int1d(start_index)].[element.field] = [neighbour_init.TradeQueueRegion][part].[element.field]
              end
            end)];        
          start_index = start_index + 1
          [neighbour_init.TradeQueueRegion][part].neighbour_part_space._valid = false
      end
    end 
    --Once we have the contiguous part of the array, we need to add neighbour_init.padding back into each cell
    for x=0, [variables.config][0].neighbour_config.x_cells do
      for y=0, [variables.config][0].neighbour_config.y_cells do
        for z=0, [variables.config][0].neighbour_config.z_cells do
           var cell_i3d = int3d({x,y,z})
          for index = 0, 50 do
            [neighbour_init.padded_particle_array][start_index + (z + y*[variables.config][0].neighbour_config.z_cells +
            x * [variables.config][0].neighbour_config.y_cells*[variables.config][0].neighbour_config.z_cells)*50 + index].neighbour_part_space._valid = false
            [neighbour_init.padded_particle_array][start_index + (z + y*[variables.config][0].neighbour_config.z_cells +
            x * [variables.config][0].neighbour_config.y_cells*[variables.config][0].neighbour_config.z_cells)*50 + index].neighbour_part_space.cell_id = cell_i3d
          end
        end
      end
    end
    --Now the padding should exist for every cell again, so we repartition. The cell mesh never changes 
    var x_cells = [variables.config][0].neighbour_config.x_cells
    var y_cells = [variables.config][0].neighbour_config.y_cells
    var z_cells = [variables.config][0].neighbour_config.z_cells
    var space_parameter = ispace(int3d, {x_cells, y_cells, z_cells}, {0,0,0})
    --Delete the old partition to reduce memory leaks
    __delete([neighbour_init.cell_partition])
    var raw_lp1 = __raw(partition([neighbour_init.padded_particle_array].neighbour_part_space.cell_id, space_parameter))
    var [neighbour_init.cell_partition] = __import_partition(disjoint, [neighbour_init.padded_particle_array], space_parameter, raw_lp1)

    --Once we've repartitioned, reintroduce the contiguity in the cells
    for cell in [neighbour_init.cell_partition].colors do
      squash_cell([neighbour_init.cell_partition][cell])
    end

    regentlib.assert(false, "Not yet implemented repartitioning and saving strategy for when tradequeues are not valid")
  end
end


return update_cells_quote
end


--Initialise the data structures needed for the TradeQueue implementation of cell lists
--The assumption is that the variables.particle_array has been initialised using IO module
--or similar and is ready to start the simulation
--TODO: Ideally after this function is called variables.particle_array
function neighbour_init.initialise(variables)

local task squash_cell(cell : region(ispace(int1d), part) ) where reads(cell), writes(cell) do

  var lo = cell.ispace.bounds.lo
  var hi = cell.ispace.bounds.hi
  var last_valid = hi
  
  for x= int32(hi), int32(lo)-1, -1 do
    if(not cell[int1d(x)].neighbour_part_space._valid) then
      --If we thought this was valid but isn't, then move the valid counter backwards
      if int1d(x) == last_valid then
        last_valid = last_valid - int1d(1)
      else
      --Otherwise swap it with the last valid entry
        [part_structure:map(function(element)
          return rquote
            cell[x].[element.field] = cell[last_valid].[element.field]
          end
        end)];
        --Set the last valid we copied to invalid
        cell[last_valid].neighbour_part_space._valid = false
        var changed = false
        -- Work backwards to find the new last valid, we know that x is always valid
          for z = int32(last_valid-1), x-1, -1 do
            --Once we find then we don't continue to search - no break statement in Lua - avoiding using goto
            if(not changed) then
              if cell[int1d(z)].neighbour_part_space._valid then
                last_valid = int1d(z)
                changed = true
              end
            end
          end
      end
    end
  end

end


--If we're running with debugging then we add an assertion to check for contiguity
local contig_assertion = function()
return rquote end
end
if DEBUG ~= nil then
local task check_contiguity(cell : region(ispace(int1d), part)) where reads(cell) do
    --Keep track of whether the tradequeues succeed
    var valid = true
    --Find end of array
    var hi = int32(cell.ispace.bounds.hi)
    --Find start of array
    var lo = int32(cell.ispace.bounds.lo)
    var first_invalid = lo
    --Find the first index containing no valid particle in the cell
    __forbid(__vectorize)
    for part = lo, hi+1 do
      if cell[int1d(part)].neighbour_part_space._valid then
        regentlib.assert(first_invalid == part, "The particles in this cell are not contiguous which is unexpected.")
        first_invalid = part + 1
      end
    end
end

contig_assertion = function()
local rval = rquote
    for cell in [neighbour_init.cell_partition].colors do
      check_contiguity([neighbour_init.cell_partition][cell])
    end
  end
  return rval
end
end


local initialisation_quote = rquote
  --Initialise the particle to cell mapping as for normal cell lists. We need to know how much to pad the array by initially
  initialise_cells(variables.config, variables.particle_array)
  particles_to_cell_launcher(variables.particle_array, variables.config)
  var num_cells = [variables.config][0].neighbour_config.x_cells *
                  [variables.config][0].neighbour_config.y_cells *
                  [variables.config][0].neighbour_config.z_cells
  --Compute the size of the particle region
  var extra_particles = neighbour_init.padding_per_cell * num_cells
  var tot_parts = extra_particles + ([variables.particle_array].ispace.bounds.hi - [variables.particle_array].ispace.bounds.lo + 1)
  var [neighbour_init.padded_particle_array] = region(ispace(int1d, tot_parts), part)
  --TODO: Initialise particle values to zero.
  --Copy the particle data into the new particle array
  for i in [variables.particle_array].ispace do
    --Generate copies here.
    [part_structure:map(function(element)
      return rquote
        [neighbour_init.padded_particle_array][i].[element.field] = [variables.particle_array][i].[element.field]
      end
    end)];
    [neighbour_init.padded_particle_array][i].neighbour_part_space._valid = true
  end
  var start_index = [variables.particle_array].ispace.bounds.hi + 1
  for x=0, [variables.config][0].neighbour_config.x_cells do
    for y=0, [variables.config][0].neighbour_config.y_cells do
      for z=0, [variables.config][0].neighbour_config.z_cells do
         var cell_i3d = int3d({x,y,z})
        for index = 0, 50 do
          [neighbour_init.padded_particle_array][start_index + (z + y*[variables.config][0].neighbour_config.z_cells +
          x * [variables.config][0].neighbour_config.y_cells*[variables.config][0].neighbour_config.z_cells)*50 + index].neighbour_part_space._valid = false
          [neighbour_init.padded_particle_array][start_index + (z + y*[variables.config][0].neighbour_config.z_cells +
          x * [variables.config][0].neighbour_config.y_cells*[variables.config][0].neighbour_config.z_cells)*50 + index].neighbour_part_space.cell_id = cell_i3d
        end
      end
    end
  end
  --We clear out the originally allocated memory because we don't really want that to exist.
  __delete([variables.particle_array])

  var n_cells = [variables.config][0].neighbour_config.x_cells * [variables.config][0].neighbour_config.y_cells * [variables.config][0].neighbour_config.z_cells
  var x_cells = [variables.config][0].neighbour_config.x_cells
  var y_cells = [variables.config][0].neighbour_config.y_cells
  var z_cells = [variables.config][0].neighbour_config.z_cells
  var space_parameter = ispace(int3d, {x_cells, y_cells, z_cells}, {0,0,0})
  var raw_lp1 = __raw(partition([neighbour_init.padded_particle_array].neighbour_part_space.cell_id, space_parameter))
  var [neighbour_init.cell_partition] = __import_partition(disjoint, [neighbour_init.padded_particle_array], space_parameter, raw_lp1)
  --TODO: TradeQueues are currently 1D partition, ideally we need to fix this.
  --https://github.com/stfc/RegentParticleDSL/issues/55
  var [neighbour_init.TradeQueueRegion] = region(ispace(int1d, neighbour_init.tradequeue_size * n_cells ), part)
  --FIXME: Initialise tradequeue values to 0 instead of this
  for i in [neighbour_init.TradeQueueRegion].ispace do
    [part_structure:map(function(element)
      return rquote
        [neighbour_init.TradeQueueRegion][i].[element.field] = [neighbour_init.padded_particle_array][0].[element.field]
      end
    end)];
    --Set invalid
    [neighbour_init.TradeQueueRegion][i].neighbour_part_space._valid = false 
  end
  var [neighbour_init.TradeQueues] = partition(equal, [neighbour_init.TradeQueueRegion], ispace(int1d, n_cells))

  --Compact the cells to make them contiguous
  for cell in [neighbour_init.cell_partition].colors do
    squash_cell([neighbour_init.cell_partition][cell])
  end
 
  for cell in [neighbour_init.cell_partition].colors do
    for i in [neighbour_init.cell_partition][cell].ispace do
      if [neighbour_init.cell_partition][cell][i].neighbour_part_space._valid and [neighbour_init.cell_partition][cell][i].core_part_space.id == int1d(8) then
        format.println("PART 8 IN CELL {} {} {}",cell.x, cell.y, cell.z)
      end
    end
  end
  --Check contiguity if debugging enabled
  [contig_assertion()]; 
end
--TODO: This might not actually work, we need to test it out.
--variables.particle_array = neighbour_init.padded_particle_array
return initialisation_quote
  
end


return neighbour_init

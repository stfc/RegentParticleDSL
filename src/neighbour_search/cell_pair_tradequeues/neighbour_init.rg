import "regent"

require("defaults")
require("src/neighbour_search/cell_pair_tradequeues/cell")

local string_to_field_path = require("src/utils/string_to_fieldpath")
local recursive_fields = require("src/utils/recursive_fields")
local neighbour_init = {}
local c = regentlib.c

get_args = require("src/utils/read_args")

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

--Read in some parameters from the command line. Work in progress
local function check_command_line()
  local debug = get_args.get_optional_arg("-hl:debug")
  if debug == 1 then
    DEBUG = true
  end
  local ppc = get_args.get_optional_arg("-hl:padding")
  if ppc ~= nil then
    neighbour_init.padding_per_cell = tonumber(ppc)
  end
  local tqsize = get_args.get_optional_arg("-hl:tqsize")
  if tqsize ~= nil then
    neighbour_init.tradequeue_size = tonumber(tqsize)
  end
  print("padding per cell is:", neighbour_init.padding_per_cell)
  print("Tradequeue size is:", neighbour_init.tradequeue_size) 
end
check_command_line()

--This function updates the cells to reflect any motion that occurs in the system. We repartition only as required, but never change
--the number of cells at this point due to causing issues with various assumptions in the system at the moment.
function neighbour_init.update_cells(variables)

--__demand(__leaf)
local task compute_new_dests(particles : region(ispace(int1d), part), config : region(ispace(int1d), config_type)) where
  reads(particles, config), writes(particles.neighbour_part_space.cell_id, particles.core_part_space) do
  for particle in particles do
    --Ignore non-valid particles
    if (particles[particle].neighbour_part_space._valid) then
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
    var tradequeue_lo = int32(tradequeue.ispace.bounds.lo)
    var tradequeue_hi = int32(tradequeue.ispace.bounds.hi)
    regentlib.assert( (tradequeue_hi-tradequeue_lo+1) == int32(tradequeue.volume), "volume != indexing size")
    var tradequeue_added = 0
    --Loop over the particles
    for part in parts.ispace do
      if(parts[part].neighbour_part_space._valid) then
        --If we're moving this particle to the trade queue
        if parts[part].neighbour_part_space.cell_id == neighbour then 
          --Check there's space in the tradequeue
          if(tradequeue_added < tradequeue.volume) then
           [part_structure:map(function(element)
              return rquote
                tradequeue[int1d(tradequeue_lo + tradequeue_added)].[element.field] = parts[part].[element.field]
              end
            end)];
            tradequeue_added = tradequeue_added + 1
            --No contiguity so we just mark this as invalid
            parts[part].neighbour_part_space._valid = false
          else
            valid = false
          end
        end
      end
    end
    --Return if successful
--    if not valid then
--      format.println("PUSH: Cell {} {} {} not valid", cell_id.x, cell_id.y, cell_id.z)
--    end
    return valid
end

--If we're running with debugging then we add an assertion to check for contiguity
local contig_assertion = function()
return rquote end
end
if DEBUG ~= nil and DEBUG then
contig_assertion = function(first_invalid, part, parts, cell)
local rval = rquote
  end
  return rval
end
end 

--__demand(__leaf)
local task tradequeue_pull(parts : region(ispace(int1d), part), tradequeue : region(ispace(int1d), part), cell : int3d) : bool where
    reads(parts, tradequeue), writes(tradequeue, parts) do

    --Keep track of whether the tradequeues succeed
    var valid = true
    --Pull some tradequeue values and check that the volume is exactly equal to the difference in index between hi and lo
    var tradequeue_lo = int32(tradequeue.ispace.bounds.lo)
    var tradequeue_hi = int32(tradequeue.ispace.bounds.hi)
    regentlib.assert( (tradequeue_hi-tradequeue_lo+1) == int32(tradequeue.volume), "volume != indexing size")
    var tradequeue_pulled = 0
    var tradequeue_remain : bool = tradequeue[tradequeue_lo].neighbour_part_space._valid
    --Loop through the cell to find empty spaces
    for part in parts.ispace do
      --If we still need to copy and we find for an empty slot
      if tradequeue_remain and not (parts[part].neighbour_part_space._valid) then
        --Copy this particle into this cell
           [part_structure:map(function(element)
              return rquote
                parts[part].[element.field] = tradequeue[int1d(tradequeue_lo + tradequeue_pulled)].[element.field]
              end
            end)];
        --Mark this tradequeue slot as empty
        tradequeue[int1d(tradequeue_lo + tradequeue_pulled)].neighbour_part_space._valid = false
        --Move to next element
        tradequeue_pulled = tradequeue_pulled + 1
        --Check if we still have particles to pull
        tradequeue_remain = tradequeue[tradequeue_lo + tradequeue_pulled].neighbour_part_space._valid
      end

    end
    --If we reach the end and tradequeue_remain is true, then we didn't have enough slots
    valid = not tradequeue_remain
--    if cell == int3d({0,0,0}) then
--      format.println("Cell 0,0,0 pulled {}", tradequeue_pulled)
--    end
--    if not valid then 
--      format.println("PULL: Cell {} {} {} not valid, pulled {}", cell.x, cell.y, cell.z, tradequeue_pulled)
--      var remaining = 0
--      for index in tradequeue.ispace do
--        if tradequeue[index].neighbour_part_space._valid then
--          remaining = remaining + 1
--        end
--      end
--      var padding = 0
--      for part in parts.ispace do
--        if not parts[part].neighbour_part_space._valid then
--          padding = padding + 1
--        end
--      end
--      format.println("Should have pulled {}, pulled {}, padding {}", remaining, tradequeue_pulled, padding)
--    end
    return valid
end

local assert_correct_cells = function()
  return rquote
  end
end
if DEBUG ~= nil and DEBUG then
  assert_correct_cells = function()
    local rval = rquote
      for cell in neighbour_init.cell_partition.colors do
        for part in neighbour_init.cell_partition[cell] do
          if neighbour_init.cell_partition[cell][part].neighbour_part_space._valid then
            regentlib.assert(neighbour_init.cell_partition[cell][part].neighbour_part_space.cell_id == int3d(cell), "particle found in wrong cell")
          end
        end
      end    
    end
    return rval
  end
end

local assert_padding = function()
  return rquote
  end
end
if DEBUG ~= nil and DEBUG then
  assert_padding = function()
    local rval = rquote
      for cell in neighbour_init.cell_partition.colors do
        var padding_count = 0
        var size = 0
        for part in neighbour_init.cell_partition[cell] do
          if not neighbour_init.cell_partition[cell][part].neighbour_part_space._valid then
            padding_count = padding_count + 1
          end
          size = size + 1
        end
--        if padding_count ~= neighbour_init.padding_per_cell then
--          format.println("Cell {} {} {} has {} padding and {} size", cell.x, cell.y, cell.z, padding_count, size)
--        end
        --if cell == int3d({0,0,0}) then
        --  format.println("Cell (0,0,0) has {} padding", padding_count)
        --end
        --regentlib.assert(padding_count == neighbour_init.padding_per_cell, "No padding") 
        regentlib.assert(padding_count ~= 0, "No padding") 
      end    
    end
  return rval
  end
end

local update_cells_quote = rquote

  var valid : bool = true;
  [assert_correct_cells()];
  [assert_padding()];
  c.legion_runtime_issue_execution_fence(__runtime(), __context())
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
  c.legion_runtime_issue_execution_fence(__runtime(), __context())
  --Trade queues attempted
  if(valid) then
    [assert_correct_cells()];
  else
  --TODO: This may need serious optimisation for occasions where it is called regularly, however we hope this never happens "too often".
  --If this is a major issue for your code, please let us know on the issue:
  --https://github.com/stfc/RegentParticleDSL/issues/56
  var tradequeue_remain = false
  var tradequeue_first : int1d
  --Loop through the tradequeue to find the index of the first particle
  for i in neighbour_init.TradeQueueRegion.ispace do
    if neighbour_init.TradeQueueRegion[i].neighbour_part_space._valid and not (tradequeue_remain) then
       tradequeue_first = i
       tradequeue_remain = true
    end
  end
  --Loop through the particles to find spaces
  for part in neighbour_init.padded_particle_array.ispace do
     --If we stil need to move things and find an empty slot
    if tradequeue_remain and not [neighbour_init.padded_particle_array][part].neighbour_part_space._valid then
      --Copy particle
           [part_structure:map(function(element)
              return rquote
                [neighbour_init.padded_particle_array][part].[element.field] = neighbour_init.TradeQueueRegion[tradequeue_first].[element.field]
              end
            end)];
           neighbour_init.TradeQueueRegion[tradequeue_first].neighbour_part_space._valid = false
           --Find another lost particle
           tradequeue_remain = false
           for i= int32(tradequeue_first), int32(neighbour_init.TradeQueueRegion.ispace.bounds.hi+1) do
             if neighbour_init.TradeQueueRegion[i].neighbour_part_space._valid and not (tradequeue_remain) then
               tradequeue_first = i
               tradequeue_remain = true
             end
           end
    end
  end

  --Now we need to repad the arrays.
  var x_cell = 0
  var y_cell = 0
  var z_cell = 0
  var counter = 0
  --Loop through the particles to find invalid entries
  for part in neighbour_init.padded_particle_array.ispace do
    --Find a space
    if not [neighbour_init.padded_particle_array][part].neighbour_part_space._valid then
      --Check we're not somehow overpadding the array.
      regentlib.assert(x_cell < [variables.config][0].neighbour_config.x_cells, "Overflowed x_cells")
      --Set it to the current cell
      var cell = int3d({x_cell, y_cell, z_cell});
      [neighbour_init.padded_particle_array][part].neighbour_part_space.cell_id = cell
      counter = counter + 1
      if counter == neighbour_init.padding_per_cell then
        counter = 0
--        format.println("Added {} padding to cell {} {} {}", neighbour_init.padding_per_cell, cell.x, cell.y, cell.z);
        z_cell = z_cell + 1
        if z_cell == [variables.config][0].neighbour_config.z_cells then
          z_cell = 0
          y_cell = y_cell + 1
          if y_cell == [variables.config][0].neighbour_config.y_cells then
            y_cell = 0
            x_cell = x_cell + 1
          end
        end
      end 
    end
  end

    --Now the padding should exist for every cell again, so we repartition. The cell mesh never changes 
    --FIXME: Do we need to delete the old partition to reduce memory leaks? At the moment doing so results in a Legion error 482
--    __delete([neighbour_init.cell_partition])
    --c.legion_runtime_issue_execution_fence(__runtime(), __context())
    --format.println("Creating new partition")
    --var n_cells = [variables.config][0].neighbour_config.x_cells * [variables.config][0].neighbour_config.y_cells * [variables.config][0].neighbour_config.z_cells
    --var x_cells = [variables.config][0].neighbour_config.x_cells
    --var y_cells = [variables.config][0].neighbour_config.y_cells
    --var z_cells = [variables.config][0].neighbour_config.z_cells
    --var space_parameter = ispace(int3d, {x_cells, y_cells, z_cells}, {0,0,0})
    --var raw_lp1 = __raw(partition([neighbour_init.padded_particle_array].neighbour_part_space.cell_id, space_parameter));
    --[neighbour_init.cell_partition] = __import_partition(disjoint, [neighbour_init.padded_particle_array], space_parameter, raw_lp1);
    --c.legion_runtime_issue_execution_fence(__runtime(), __context());
    ----Check partition happened correctly
    --[assert_correct_cells()];
    --[assert_padding()];
    --valid = true
    regentlib.assert(false, "Not yet implemented repartitioning and saving strategy for when tradequeues are not valid. Use command line flags -hl:padding (default 50) and -hl:tqsize (default 10 to increase the extra space to allow particle congestion")
  end
end


return update_cells_quote
end


--Initialise the data structures needed for the TradeQueue implementation of cell lists
--The assumption is that the variables.particle_array has been initialised using IO module
--or similar and is ready to start the simulation
--TODO: Ideally after this function is called variables.particle_array
function neighbour_init.initialise(variables)


--If we're running with debugging then we add an assertion to check for contiguity
local contig_assertion = function()
return rquote end
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
        for index = 0, neighbour_init.padding_per_cell do
          [neighbour_init.padded_particle_array][start_index + (z + y*[variables.config][0].neighbour_config.z_cells +
          x * [variables.config][0].neighbour_config.y_cells*[variables.config][0].neighbour_config.z_cells)*neighbour_init.padding_per_cell + index].neighbour_part_space._valid = false
          [neighbour_init.padded_particle_array][start_index + (z + y*[variables.config][0].neighbour_config.z_cells +
          x * [variables.config][0].neighbour_config.y_cells*[variables.config][0].neighbour_config.z_cells)*neighbour_init.padding_per_cell + index].neighbour_part_space.cell_id = cell_i3d
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
  var raw_lp1 = __raw(partition([neighbour_init.padded_particle_array].neighbour_part_space.cell_id, space_parameter));
  var [neighbour_init.cell_partition] = __import_partition(disjoint, [neighbour_init.padded_particle_array], space_parameter, raw_lp1)
  --TODO: TradeQueues are currently 1D partition, ideally we need to fix this.
  --https://github.com/stfc/RegentParticleDSL/issues/55
  var [neighbour_init.TradeQueueRegion] = region(ispace(int1d, neighbour_init.tradequeue_size * n_cells ), part)
  --FIXME: Initialise tradequeue values to 0 instead of this - requires zero_parts functionality
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

end
return initialisation_quote
  
end


return neighbour_init
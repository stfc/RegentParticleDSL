include "regent"

require("defaults")
require("src/neighbour_search/cell_pair_tradequeues/cell")

local string_to_field_path = require("src/utils/string_to_fieldpath")
local recursive_fields = require("src/utils/recursive_fields")
local neighbour_init = {}

neighbour_init.padded_particle_array = regentlib.newsymbol("padded_particle_array")
neighbour_init.TradeQueues = regentlib.newsymbol("TradeQueues")

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

__demand(__leaf)
local task compute_new_dests(parts : region(ispace(int1d), part), config : region(ispace(int1d), config)) where
  reads(parts, config), writes(parts.cell_id, parts.core_part_space) do
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

__demand(__leaf)
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
    --Inside terra so C style loop, x=hi, x > lo, x--
    for part=hi, lo-1, -1 do
      if(parts[int1d(part)]._valid) then
        --If we're moving this particle to the trade queue
        if parts[int1d(part)].neighbour_part_space.cell_id = neighbour then
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
contig_assertion = function(first_invalid, part)
local rval = rquote
  regentlang.assert(first_invalid == part-1, "The particles in this cell are not contiguous which is unexpected.")
  end
  return rval
end
end 
__demand(__leaf)
local task tradequeue_pull(parts : region(ispace(int1d), part), tradequeue : region(ispace(int1d), part)) : bool where
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
        [contig_assertion(first_invalid, part)];
        first_invalid = part + 1
      end 
    end
    for part in tradequeue.ispace do
      if valid and tradequeue[part].neighbour_part_space._valid then
        if first_invalid <= hi then
          --Copy the particle data into the first invalid index
            [part_structure:map(function(element)
              return rquote
                parts[first_invalid].[element.field] = tradequeue[part].[element.field]
              end
            end)];
          first_invalid = first_invalid + int1d(1)
          tradequeue[part]._valid = false
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
      for cell in variables.cell_partitions.colors do
        for part in variables.cell_partition[cell] do
          regentlib.assert(variables.cell_partition[cell][part].neighbour_part_space.cell_id == cell, "Cell found in wrong cell")
        end
      end    
    end
    return rval
  end
end

local update_cells_quote = rquote

  var valid : bool = true
  __demand(__index_launch)
  for cell in variables.cell_partition.colors do
    compute_new_dests([variables.cell_partition][cell])
  end
  for direction = 0, 26 do
    if(valid) then
      --Push and pull tradequeues
      --Compute the neighbour in the relevant direction
      var new_x = cell.x + directions[direction][0]
      var new_y = cell.y + directions[direction][1]
      var new_z = cell.z + directions[direction][2]
      if new_x < 0 then new_x = new_x + [variables.config][0].neighbour_config.x_cells end
      if new_y < 0 then new_y = new_y + [variables.config][0].neighbour_config.y_cells end
      if new_z < 0 then new_z = new_z + [variables.config][0].neighbour_config.z_cells end
      if new_x > [variables.config][0].neighbour_config.x_cells then new_x = new_x - [variables.config][0].neighbour_config.x_cells end
      if new_y > [variables.config][0].neighbour_config.y_cells then new_y = new_y - [variables.config][0].neighbour_config.y_cells end
      if new_z > [variables.config][0].neighbour_config.z_cells then new_z = new_z - [variables.config][0].neighbour_config.z_cells end
      var neighbour : int3d = int3d({new_x, new_y, new_z})
      --First spawn all the push tasks
      for cell in variables.cell_partitions.colors do
        valid = valid and tradequeue_push([variables.cell_partition][cell], neighbour_init.TradeQueues[cell], cell, neighbour)
      end
      --Once all the push tasks are spawned, spawn the pull tasks
      for cell in variables.cell_partitions.colors do
        valid = valid and tradequeue_pull([variables.cell_partition][neighbour], neighbour_init.TradeQueues[cell])
      end
    end
  end
  --Trade queues attempted
  if(valid) then
    [assert_correct_cells()];
  else
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

local initalisation_quote = rquote
  --Initialise the particle to cell mapping as for normal cell lists. We need to know how much to pad the array by initially
  initialise_cells(variables.config, variables.particle_array)
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
  for x=0, x_cells do
    for y=0, y_cells do
      for z=0, z_cells do
        [neighbour_init.padded_particle_array][start_index + z + y*z_cells + x * y_cells*z_cells]._valid = false
      end
    end
  end
  --We clear out the originally allocated memory because we don't really want that to exist.
  __delete([variables.particle_array])

  --TODO: Partition the particles
  --TODO: Initialise the trade queues
  
end
--TODO: This might not actually work, we need to test it out.
variables.particle_array = neighbour_init.padded_particle_array
return initialisation_quote
  
end


return neighbour_init

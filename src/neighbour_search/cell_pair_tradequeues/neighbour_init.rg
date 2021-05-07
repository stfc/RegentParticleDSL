import "regent"

require("src/neighbour_search/cell_pair_tradequeues/cell")
require("src/particles/init_part")

local string_to_field_path = require("src/utils/string_to_fieldpath")
local recursive_fields = require("src/utils/recursive_fields")
local neighbour_init = {}
local c = regentlib.c

get_args = require("src/utils/read_args")

neighbour_init.padded_particle_array = regentlib.newsymbol("padded_particle_array")
---Maybe no longer need a TradeQueueRegion separately.
--neighbour_init.TradeQueueRegion = regentlib.newsymbol("TradeQueueRegion")
--Create 26 tradequeues similar to Soleil-X
neighbour_init.TradeQueues = terralib.newlist()
for i=1, 26 do
 --neighbour_init.TradeQueues:insert( regentlib.newsymbol(region(ispace(int1d), part) ) )
 neighbour_init.TradeQueues:insert( regentlib.newsymbol( ) )
end

--Create symbols for the src/dest partitions
neighbour_init.TradeQueues_bySrc = terralib.newlist()
for i=1, 26 do
  neighbour_init.TradeQueues_bySrc:insert( regentlib.newsymbol() )
end

neighbour_init.TradeQueues_byDest = terralib.newlist()
for i=1, 26 do
  neighbour_init.TradeQueues_byDest:insert( regentlib.newsymbol() )
end

neighbour_init.cell_partition = regentlib.newsymbol("padded_cell_partition")

local DEBUG = true

--Use terralib list for this isntead of lua table
local directions = terralib.newlist({
    rexpr int3d({-1, -1, -1}) end,
    rexpr int3d({-1, -1,  0}) end,
    rexpr int3d({-1, -1,  1}) end,
    rexpr int3d({-1,  0, -1}) end,
    rexpr int3d({-1,  0,  0}) end,
    rexpr int3d({-1,  0,  1}) end,
    rexpr int3d({-1,  1, -1}) end,
    rexpr int3d({-1,  1,  0}) end,
    rexpr int3d({-1,  1,  1}) end,
    rexpr int3d({ 0, -1, -1}) end,
    rexpr int3d({ 0, -1,  0}) end,
    rexpr int3d({ 0, -1,  1}) end,
    rexpr int3d({ 0,  0, -1}) end,
    rexpr int3d({ 0,  0,  1}) end,
    rexpr int3d({ 0,  1, -1}) end,
    rexpr int3d({ 0,  1,  0}) end,
    rexpr int3d({ 0,  1,  1}) end,
    rexpr int3d({ 1, -1, -1}) end,
    rexpr int3d({ 1, -1,  0}) end,
    rexpr int3d({ 1, -1,  1}) end,
    rexpr int3d({ 1,  0, -1}) end,
    rexpr int3d({ 1,  0,  0}) end,
    rexpr int3d({ 1,  0,  1}) end,
    rexpr int3d({ 1,  1, -1}) end,
    rexpr int3d({ 1,  1,  0}) end,
    rexpr int3d({ 1,  1,  1}) end,
})

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

--Find where the particles now belong
--This task is still lightweight and could maybe combined into TradeQueue_push
local __demand(__leaf) task compute_new_dests(particles : region(ispace(int1d), part), config : region(ispace(int1d), config_type)) where
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

--Tradequeue symbols used for the generation of the push/pull tasks
local tradequeues = terralib.newlist()
for i=1, 26 do
  tradequeues:insert( regentlib.newsymbol(region(ispace(int1d), part) ) )
end

--We need to add a flatten function for the tradequeue flattening to work.
--FIXME: This needs redoing to be better/distinct
local TerraList = getmetatable(terralib.newlist())
function TerraList:flatten(res)
    res = res or terralib.newlist()
    for _, v in pairs(self) do
        if terralib.israwlist(v) then
           v:flatten(res)
        else
            res:insert(v)
        end
    end
    return res
end

--New tradequeue push implementation
local __demand(__leaf) task tradequeue_push( cell_id : int3d, 
                                             particles : region(ispace(int1d), part),
                                             config : region(ispace(int1d), config_type),
                                             [tradequeues] ) where
      reads(particles), reads(config), reads writes(particles.neighbour_part_space),
      --We write to all the fields in each of the tradeQueue partitions
      [tradequeues:map(function(queue)
        return part_structure:map(function (element)
            return regentlib.privilege(regentlib.writes, queue, element.field)
        end)
      end):flatten()]
do

    var toTransfer = int32(0)
    var count_xcells = config[0].neighbour_config.x_cells
    var count_ycells = config[0].neighbour_config.y_cells
    var count_zcells = config[0].neighbour_config.z_cells
    for part in particles.ispace do
        --We only care about valid particles
        if particles[part].neighbour_part_space._valid then
            --Check if it still belongs to us
            if particles[part].neighbour_part_space.cell_id ~= cell_id then
                toTransfer += 1
                --Find out which neighbour it moves to
                [(function() local __quotes = terralib.newlist()
                  for i = 1, 26 do
                      __quotes:insert(rquote
                        if particles[part].neighbour_part_space._transfer_dir == 0 and
                           (particles[part].neighbour_part_space.cell_id
                            == (cell_id + [directions[i]] + {count_xcells,count_ycells,count_zcells})%{count_xcells,count_ycells,count_zcells} )
                           then
                            particles[part].neighbour_part_space._transfer_dir = i
                            toTransfer -= 1
                        end
                      end)
                  end
                  return __quotes
                 end) ()];
            end
        end
    end
    if toTransfer ~= 0 then
        var s : rawstring
        s = [rawstring] (regentlib.c.malloc(512))
        for part in particles.ispace do
            if particles[part].neighbour_part_space.cell_id ~= cell_id and particles[part].neighbour_part_space._transfer_dir == 0 then
                format.println("Part should go to {} {} {}  from {} {} {} and has velocity {} {} {}", 
                                particles[part].neighbour_part_space.cell_id.x,
                                particles[part].neighbour_part_space.cell_id.y,
                                particles[part].neighbour_part_space.cell_id.z,
                                cell_id.x, cell_id.y, cell_id.z, particles[part].core_part_space.vel_x,
                                particles[part].core_part_space.vel_y, particles[part].core_part_space.vel_z)
            end
        end
        format.snprintln(s,512, "Particle moving to non-adjacent cell from cell {} {} {}", cell_id.x, cell_id.y, cell_id.z)
        regentlib.assert(toTransfer == 0, s)
        regentlib.c.free(s)
    end

   --Loop through the queues, clear them and copy the particle data. Each cell owns its own section of the queue for each direction so
   --don't need to worry about any other cell touching it.
    [(function() local __quotes = terralib.newlist()
        for i=1, 26 do
           local queue = tradequeues[i]
            __quotes:insert(rquote
                --Empty the queue from last step
                for j in queue.ispace do
                    queue[j].neighbour_part_space._valid = false
                end
               --Count how many to transfer and where they belong
                var transferred = int32(0)
                for part in particles.ispace do
                    if particles[part].neighbour_part_space._valid and particles[part].neighbour_part_space._transfer_dir == i then
                        particles[part].neighbour_part_space._transfer_pos = int1d(transferred)
                        transferred += 1
                    else
                        particles[part].neighbour_part_space._transfer_pos = int1d(0)
                    end
                end
                --Check we have enough space
                if int1d(transferred) > queue.bounds.hi - queue.bounds.lo + 1 then
                    var s : rawstring
                    s = [rawstring] (regentlib.c.malloc(512))
                    format.snprintln(s, 512, "Transferring more particles than fit in the tradequeue. Cell {} {} {}, transfers: {}", 
                                     cell_id.x, cell_id.y, cell_id.z, int32(transferred))
                    regentlib.assert(int1d(transferred) <= (queue.bounds.hi - queue.bounds.lo + 1), s)
                    regentlib.c.free(s)
                end
                --We have enough space, so move the particles into the queue
                for part in particles.ispace do
                    if particles[part].neighbour_part_space._valid and particles[part].neighbour_part_space._transfer_dir == i then
                        var pos = particles[part].neighbour_part_space._transfer_pos + queue.bounds.lo;
                        --Copy the particle data.
                        [part_structure:map(function(element)
                            return rquote
                                queue[pos].[element.field] = particles[part].[element.field]
                            end
                         end)];
                        --Mark the particle no longer valid
                        particles[part].neighbour_part_space._valid = false
                    end
                end
            end)
        end
        return __quotes
    end) ()];
    --All done now
end


--: New tradequeue pull implementation
local __demand(__leaf) task tradequeue_pull( cell_id : int3d,
                                             particles : region(ispace(int1d), part),
                                             config : region(ispace(int1d), config_type),
                                             [tradequeues] ) where
      reads(config), reads writes(particles),
      --We read to all the fields in each of the tradeQueue partitions
      [tradequeues:map(function(queue)
        return part_structure:map(function (element)
            return regentlib.privilege(regentlib.reads, queue, element.field)
        end) 
      end):flatten()]
do
    --Keep track of where to transfer elements and total number of transfers
    var transfer_bounds : int32[27]
    transfer_bounds[0] = 0
    var total_transfers : int32 = 0;
    --Count the number of transfers
    [(function() local __quotes = terralib.newlist()
        for i=1, 26 do
           local queue = tradequeues[i]
            __quotes:insert(rquote
                for q in queue.ispace do
                    if queue[q].neighbour_part_space._valid then
                        total_transfers += 1 
                    end
                end
                transfer_bounds[i] = total_transfers
            end)
        end
        return __quotes
    end) ()];
    --Number all the slots and check we have enough space
    var available_slots : int32 = 0
    for part in particles.ispace do
        if not (particles[part].neighbour_part_space._valid) then
            particles[part].neighbour_part_space._transfer_pos = int1d(available_slots)
            available_slots += 1
        end
    end
    --Check there's enough space in the particle array to accept the particles
    if available_slots < total_transfers then
        var s : rawstring
        s = [rawstring] (regentlib.c.malloc(512))
        format.snprintln(s,512, "Not enough space in cell {} {} {} to accept incoming particles. {} space and {} required", cell_id.x, cell_id.y, cell_id.z,
                                available_slots, total_transfers)
        regentlib.assert(available_slots >= total_transfers, s)
        regentlib.c.free(s)
    end
    --Ok, we're good to go!
    [(function() local __quotes = terralib.newlist()
        for i=1, 26 do
            local queue = tradequeues[i]
            __quotes:insert(rquote
                --Get the bounds for this direction
                var lo = int1d(transfer_bounds[i-1])
                var hi = int1d(transfer_bounds[i])
                for part in particles.ispace do
                    if not(particles[part].neighbour_part_space._valid) then
                        var offset = particles[part].neighbour_part_space._transfer_pos
                        if offset >= lo and offset < hi then
                            var pos : int1d = offset - lo + queue.bounds.lo;
                            [part_structure:map(function(element)
                               return rquote
                                 particles[part].[element.field] = queue[pos].[element.field]
                               end
                             end)];                            
                            particles[part].neighbour_part_space._transfer_dir = 0
                            particles[part].neighbour_part_space._transfer_pos = 0
                        end
                    end
                end            
            end)
        end
        return __quotes
    end) ()];
end

--Partitioning code for the tradequeues to correctly create the src/dest tradequeues.
local __demand(__inline) task partition_tradequeue_by_cells( tradequeue : region(ispace(int1d), part),
                                                             cell_space : ispace(int3d),
                                                             offset : int3d)
    var count = tradequeue.bounds.hi + 1
    var count_xcells = cell_space.bounds.hi.x + 1
    var count_ycells = cell_space.bounds.hi.y + 1
    var count_zcells = cell_space.bounds.hi.z + 1
    var n_cells = count_xcells * count_ycells * count_zcells

    --Use legion's coloring option to create this partition
    var coloring = regentlib.c.legion_domain_point_coloring_create()
    for cell in cell_space do
      --Ensure we wrap the cell indices (we add # cells to each dimension to handle negative values)
      var color : int3d = (cell - offset + {count_xcells, count_ycells, count_zcells}) % {count_xcells, count_ycells, count_zcells}
      --Indices for this color are number of elements per cell (count/n_cells) * (3D to 1D conversion)
      var oneD_color : int1d = (color.x*count_ycells*count_zcells) + (color.y*count_zcells) + color.z;
      var rect = rect1d{
        lo = (count/n_cells)*(oneD_color),
        hi = (count/n_cells)*(oneD_color+1) - 1
      }
      regentlib.c.legion_domain_point_coloring_color_domain(coloring, cell, rect)
    end
    --Create the partition from the coloring. If offset = 0 this is essentially a controlled equal partition
    var p = partition(disjoint, tradequeue, coloring, cell_space)
    --Clean up the runtime
    regentlib.c.legion_domain_point_coloring_destroy(coloring)
    return p
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


local function generate_range_as_terralist(lo, hi)

    local range = terralib.newlist()
    for i=lo, hi do
        range:insert(i)
    end
    return range
end

--local function generate_tradequeue_src_rexpr(cell)
--    local range = terralib.newlist()
--    for i=1, 26 do
--        range:insert(i)
--    end
--    return range:map(function(i) return rexpr
--            [neighbour_init.TradeQueues_bySrc[i]][cell]
--        end
--    end)
--end
--
--local function generate_tradequeue_dest_rexpr(cell)
--    local range = terralib.newlist()
--    for i=1, 26 do
--        range:insert(i)
--    end
--    return range:map(function(i) return rexpr
--            [neighbour_init.TradeQueues_byDest[i]][cell]
--        end
--    end);
--end

local function get_timing_quotes()
  local start_timing_quote = rquote

  end
  local end_timing_quote = rquote

  end
  if DSL_settings.TIMING then
    local starttime = regentlib.newsymbol()
    local endtime = regentlib.newsymbol()
    start_timing_quote = rquote
      var [starttime] = c.legion_get_current_time_in_micros()
    end
    end_timing_quote = rquote
      c.legion_runtime_issue_execution_fence(__runtime(), __context())
      var [endtime] = c.legion_get_current_time_in_micros();
      [variables.config][0].timing_config.neighbour_search_time = [variables.config][0].timing_config.neighbour_search_time + ([endtime] - [starttime])
    end
  end
  return start_timing_quote, end_timing_quote
end

--This function updates the cells to reflect any motion that occurs in the system. We repartition only as required, but never change
--the number of cells at this point due to causing issues with various assumptions in the system at the moment.
function neighbour_init.update_cells(variables)
local start_timing_quote, end_timing_quote = get_timing_quotes()
local update_cells_quote = rquote
    [start_timing_quote];
    [assert_correct_cells()];
    for cell in [neighbour_init.cell_partition].colors do
        compute_new_dests( [neighbour_init.cell_partition][cell], [variables.config]);
    end
    for cell in [neighbour_init.cell_partition].colors do
        tradequeue_push(cell, [neighbour_init.cell_partition][cell], [variables.config],
                        [generate_range_as_terralist(1, 26):map(function(i) return rexpr
                            [neighbour_init.TradeQueues_bySrc[i]][cell]
                            end
                        end)] )
    end
    for cell in [neighbour_init.cell_partition].colors do
        tradequeue_pull(cell, [neighbour_init.cell_partition][cell], [variables.config],
                    [generate_range_as_terralist(1, 26):map(function(i) return rexpr
                        [neighbour_init.TradeQueues_byDest[i]][cell]
                        end
                    end)] )
    end
    c.legion_runtime_issue_execution_fence(__runtime(), __context());
    [assert_correct_cells()];
    [end_timing_quote];
    -- TODO: If you need repartitioning for your testcase and tradequeues are not sufficient, we need a new implementation of neighbour search
    -- If this is something you need, let us know on this issue:
    -- https://github.com/stfc/RegentParticleDSL/issues/56
end

return update_cells_quote
end

--Initialise the data structures needed for the TradeQueue implementation of cell lists
--The assumption is that the variables.particle_array has been initialised using IO module
--or similar and is ready to start the simulation
--Ideally after this function is called variables.particle_array
function neighbour_init.initialise(variables)

local initialisation_quote = rquote                                                                                                                                          --Initialise the particle to cell mapping as for normal cell lists. We need to know how much to pad the array by initially
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
  __delete([variables.particle_array]);

--Init Tradequeues
  --The wrapper here is an inline way to create this for all 26 directions inside Regent code, since we're using
  -- symbols to generate this. We could do this with neighbour_init.TradeQueues:map(...) but we may in the future
  --vary tradequeue size based on direction.
  [(function() local __quotes = terralib.newlist() 
     for i = 1,26 do
     __quotes:insert(rquote
       --We could vary tradequeue size based upon direction, but for now we don't care.
       var TradeQueue_indexSpace = ispace( int1d, num_cells * neighbour_init.tradequeue_size )
       var [neighbour_init.TradeQueues[i]] = region(TradeQueue_indexSpace, part);
         [generate_zero_part_quote( neighbour_init.TradeQueues[i] )];
     end)
   end
   return __quotes
   end) ()];

  --Now we have to partition the tradequeues by source and by destination.
  --For the source we just divide each of the 26 "direction's" tradequeues into one per source cell. We use regentlib.coloring
  --code to do this so we know the size/position accurately.
  --For the destination, we do the same, however we use an offset (based upon the direction) to shift each cells tradequeue appropriately.
  var x_cells = [variables.config][0].neighbour_config.x_cells
  var y_cells = [variables.config][0].neighbour_config.y_cells
  var z_cells = [variables.config][0].neighbour_config.z_cells
  var cell_space = ispace(int3d, {x_cells, y_cells, z_cells});
  [(function() local __quotes = terralib.newlist()
    for i = 1, 26 do
        __quotes:insert(rquote
            var [neighbour_init.TradeQueues_bySrc[i]] = partition_tradequeue_by_cells( [neighbour_init.TradeQueues[i]],
                                                                                       cell_space, int3d({0,0,0}) );
            var [neighbour_init.TradeQueues_byDest[i]] = partition_tradequeue_by_cells( [neighbour_init.TradeQueues[i]],
                                                                                        cell_space, [directions[i]] );
        end)
    end
    return __quotes
   end) ()];

   --Init the cell partition
    var space_parameter = ispace(int3d, {x_cells, y_cells, z_cells}, {0,0,0})
    var raw_lp1 = __raw(partition([neighbour_init.padded_particle_array].neighbour_part_space.cell_id, space_parameter));
    var [neighbour_init.cell_partition] = __import_partition(disjoint, [neighbour_init.padded_particle_array], space_parameter, raw_lp1)
end

return initialisation_quote
end

return neighbour_init

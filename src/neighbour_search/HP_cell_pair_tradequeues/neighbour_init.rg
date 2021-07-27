import "regent"

require("src/neighbour_search/HP_cell_pair_tradequeues/cell")
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

neighbour_init.TradeQueues_bysubSrc = terralib.newlist()
for i=1, 26 do
    neighbour_init.TradeQueues_bysubSrc:insert(regentlib.newsymbol() )
end

neighbour_init.TradeQueues_bysubDest = terralib.newlist()
for i=1, 26 do
    neighbour_init.TradeQueues_bysubDest:insert(regentlib.newsymbol() )
end

neighbour_init.cell_partition = regentlib.newsymbol("padded_cell_partition")
neighbour_init.supercell_partition = regentlib.newsymbol("supercell_partition")
neighbour_init.halo_partition = regentlib.newsymbol("halo_partition")
neighbour_init.x_slices = regentlib.newsymbol("x_slice_partition")


local DEBUG = true

--Use terralib list for this isntead of lua table
--!There are 13 positive vectors to sort along.
--!These depend on the shape of the volume (since this affects the shape of the cells)
--!The 13 vectors to sort along are (in order)
--! 0 1 0, 0 0 1, 0 1 1, 0 1 -1
--! 1 0 0, 1 0 1, 1 1 0, 1 0 -1
--! 1 -1 0, 1 1 1, 1 -1 -1, 1 -1 1
--! 1 1 -1
local directions = terralib.newlist({
    rexpr int3d({ 0,  1,  0}) end,
    rexpr int3d({ 0,  0,  1}) end,
    rexpr int3d({ 0,  1,  1}) end,
    rexpr int3d({ 0,  1, -1}) end,
    rexpr int3d({ 1,  0,  0}) end,
    rexpr int3d({ 1,  0,  1}) end,
    rexpr int3d({ 1,  1,  0}) end,
    rexpr int3d({ 1,  0, -1}) end,
    rexpr int3d({ 1, -1,  0}) end,
    rexpr int3d({ 1,  1,  1}) end,
    rexpr int3d({ 1, -1, -1}) end,
    rexpr int3d({ 1, -1,  1}) end,
    rexpr int3d({ 1,  1, -1}) end,

    rexpr int3d({ 0, -1,  0}) end,
    rexpr int3d({ 0,  0, -1}) end,
    rexpr int3d({ 0, -1, -1}) end,
    rexpr int3d({ 0, -1,  1}) end,
    rexpr int3d({-1,  0,  0}) end,
    rexpr int3d({-1,  0, -1}) end,
    rexpr int3d({-1, -1,  0}) end,
    rexpr int3d({-1,  0,  1}) end,
    rexpr int3d({-1,  1,  0}) end,
    rexpr int3d({-1, -1, -1}) end,
    rexpr int3d({-1,  1,  1}) end,
    rexpr int3d({-1,  1, -1}) end,
    rexpr int3d({-1, -1,  1}) end
})
    

--local directions = terralib.newlist({
--    rexpr int3d({-1, -1, -1}) end,
--    rexpr int3d({-1, -1,  0}) end,
--    rexpr int3d({-1, -1,  1}) end,
--    rexpr int3d({-1,  0, -1}) end,
--    rexpr int3d({-1,  0,  0}) end,
--    rexpr int3d({-1,  0,  1}) end,
--    rexpr int3d({-1,  1, -1}) end,
--    rexpr int3d({-1,  1,  0}) end,
--    rexpr int3d({-1,  1,  1}) end,
--    rexpr int3d({ 0, -1, -1}) end,
--    rexpr int3d({ 0, -1,  0}) end,
--    rexpr int3d({ 0, -1,  1}) end,
--    rexpr int3d({ 0,  0, -1}) end,
--    rexpr int3d({ 0,  0,  1}) end,
--    rexpr int3d({ 0,  1, -1}) end,
--    rexpr int3d({ 0,  1,  0}) end,
--    rexpr int3d({ 0,  1,  1}) end,
--    rexpr int3d({ 1, -1, -1}) end,
--    rexpr int3d({ 1, -1,  0}) end,
--    rexpr int3d({ 1, -1,  1}) end,
--    rexpr int3d({ 1,  0, -1}) end,
--    rexpr int3d({ 1,  0,  0}) end,
--    rexpr int3d({ 1,  0,  1}) end,
--    rexpr int3d({ 1,  1, -1}) end,
--    rexpr int3d({ 1,  1,  0}) end,
--    rexpr int3d({ 1,  1,  1}) end,
--})

local function construct_part_structure()
  local part_structure = terralib.newlist()
  local field_strings = {}
  local type_table = {}
  for k, v in pairs(part.fields) do
    recursive_fields.recurse_field(v, field_strings, type_table)
  end
  for k, _ in pairs(field_strings) do
    local st = string.find(field_strings[k], "halos")
    if st == nil then
        part_structure:insert({field = string_to_field_path.get_field_path(field_strings[k])})
    end
  end
  return part_structure
end
local part_structure = construct_part_structure()

fspace sorting_ids{
    sid : int1d[13]
}
neighbour_init.sorting_id_array = regentlib.newsymbol("sorting_array")
neighbour_init.sorting_array_cell_partition = regentlib.newsymbol("sorting_array_cell_partition")
neighbour_init.sorting_array_supercell_partition = regentlib.newsymbol("sorting_array_supercell_partition")

--This realistically is simulation dependent
neighbour_init.padding_per_cell = 10
neighbour_init.tradequeue_size = neighbour_init.padding_per_cell

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
  reads(particles, config), writes(particles.neighbour_part_space.cell_id, particles.neighbour_part_space.x_cell, particles.core_part_space,
                                   particles.neighbour_part_space.supercell_id) do
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
      var x_supercell : int1d = int1d( (particles[particle].core_part_space.pos_x / config[0].neighbour_config.supercell_dim_x) )
      var y_supercell : int1d = int1d( (particles[particle].core_part_space.pos_x / config[0].neighbour_config.supercell_dim_y) )
      var z_supercell : int1d = int1d( (particles[particle].core_part_space.pos_x / config[0].neighbour_config.supercell_dim_z) )
      cell_loc = int3d( {x_supercell, y_supercell, z_supercell} )
      particles[particle].neighbour_part_space.supercell_id = cell_loc
      particles[particle].neighbour_part_space.x_cell = x_supercell
    end
  end
end

--Tradequeue symbols used for the generation of the push/pull tasks
local tradequeues = terralib.newlist()
for i=1, 26 do
  tradequeues:insert( regentlib.newsymbol(region(ispace(int1d), part) ) )
end

local all_tradequeues = terralib.newlist()
for i= 1, 26 do
    all_tradequeues:insert(regentlib.newsymbol(region(ispace(int1d), part) ) )
end

local tradequeues_partitions = terralib.newlist()
for i=1, 26 do
  tradequeues_partitions:insert(regentlib.newsymbol(partition(disjoint,all_tradequeues[i], ispace(int3d))))
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
local __demand(__leaf) task tradequeue_push( supercell_id : int3d, 
                                             particles : region(ispace(int1d), part),
                                             all_parts : region(ispace(int1d), part),
                                             subcell_partition : partition(disjoint, all_parts, ispace(int3d)),
                                             config : region(ispace(int1d), config_type),
                                             [tradequeues], [all_tradequeues], [tradequeues_partitions] ) where
      reads(particles), reads(config), reads writes(particles.neighbour_part_space),
      --We write to all the fields in each of the tradeQueue partitions
      [tradequeues:map(function(queue)
        return part_structure:map(function (element)
            return regentlib.privilege(regentlib.writes, queue, element.field)
        end)
      end):flatten()]
do

    var count_xcells = config[0].neighbour_config.x_cells
    var count_ycells = config[0].neighbour_config.y_cells
    var count_zcells = config[0].neighbour_config.z_cells
    var n_cells = count_xcells * count_ycells * count_zcells

    var x_per_super = config[0].neighbour_config.x_cells / config[0].neighbour_config.x_supercells
    var y_per_super = config[0].neighbour_config.y_cells / config[0].neighbour_config.y_supercells
    var z_per_super = config[0].neighbour_config.z_cells / config[0].neighbour_config.z_supercells

    --Compute the bounds for the subcells of this supercell
    var xlo = supercell_id.x * x_per_super
    var xhi = (supercell_id.x + 1) * x_per_super --loops noninclusive so all ok to not -1
    var ylo = supercell_id.y * y_per_super
    var yhi = (supercell_id.y + 1) * y_per_super
    var zlo = supercell_id.z * z_per_super
    var zhi = (supercell_id.z + 1) * z_per_super

--format.println("supercell {} {} {} index space:", supercell_id.x, supercell_id.y, supercell_id.z)
--for part in particles.ispace do
--    format.println("{}", part)
--end
--format.println("End of index space")

var total_transfers = int32(0);
for x=xlo, xhi do
for y=ylo, yhi do
for z=zlo, zhi do
    var cell_id : int3d = int3d({x,y,z})
--    format.println("supercell {} {} {}, subcell {} {} {}, per super {} {} {}", supercell_id.x, supercell_id.y, supercell_id.z, cell_id.x, cell_id.y, cell_id.z,
--                    x_per_super, y_per_super, z_per_super)
    --Lets do a bounds check
--    for part in subcell_partition[cell_id].ispace do
--        var found = false
--        for part2 in particles.ispace do
--            if int1d(part) == int1d(part2) then
--                found = true
--            end
--        end
--        if not found then
--            format.println("Expected to find index {} but only was not in supercell", part)
--           regentlib.assert(found, "Didn't find index")
--        end
--    end
    var toTransfer = int32(0)
    for part in subcell_partition[cell_id].ispace do
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
           local queuespace = tradequeues_partitions[i]
           local queue = tradequeues[i]
            __quotes:insert(rquote
                --Empty the queue from last step
                for j in queuespace[cell_id].ispace do
                    queue[j].neighbour_part_space._valid = false
                end
               --Count how many to transfer and where they belong
                var transferred = int32(0)
                for part in subcell_partition[cell_id].ispace do
                    if particles[part].neighbour_part_space._valid and particles[part].neighbour_part_space._transfer_dir == i then
                        particles[part].neighbour_part_space._transfer_pos = int1d(transferred)
                        transferred += 1
                    else
                        particles[part].neighbour_part_space._transfer_pos = int1d(0)
                    end
                end
                --Check we have enough space
                if int1d(transferred) > queuespace[cell_id].bounds.hi - queuespace[cell_id].bounds.lo + 1 then
                    var s : rawstring
                    s = [rawstring] (regentlib.c.malloc(512))
                    format.snprintln(s, 512, "Transferring more particles than fit in the tradequeue. Cell {} {} {}, transfers: {}, size: {}", 
                                     cell_id.x, cell_id.y, cell_id.z, int32(transferred), queuespace[cell_id].bounds.hi - queuespace[cell_id].bounds.lo)
                    regentlib.assert(int1d(transferred) <= (queuespace[cell_id].bounds.hi - queuespace[cell_id].bounds.lo + 1), s)
                    regentlib.c.free(s)
                end
                total_transfers = total_transfers + transferred
                --We have enough space, so move the particles into the queue
                for part in subcell_partition[cell_id].ispace do
                    if particles[part].neighbour_part_space._valid and particles[part].neighbour_part_space._transfer_dir == i then
                        var pos = particles[part].neighbour_part_space._transfer_pos + queuespace[cell_id].bounds.lo;
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
    --All done now for this cell
end
end
end

    --All done now
end


--: New tradequeue pull implementation
local __demand(__leaf) task tradequeue_pull( supercell_id : int3d,
                                             particles : region(ispace(int1d), part),
                                             all_parts : region(ispace(int1d), part),
                                             subcell_partition : partition(disjoint, all_parts, ispace(int3d)),
                                             config : region(ispace(int1d), config_type),
                                             [tradequeues], [all_tradequeues], [tradequeues_partitions] ) where
      reads(config), reads writes(particles),
      --We read to all the fields in each of the tradeQueue partitions
      [tradequeues:map(function(queue)
        return part_structure:map(function (element)
            return regentlib.privilege(regentlib.reads, queue, element.field)
        end) 
      end):flatten()]
do
    var count_xcells = config[0].neighbour_config.x_cells
    var count_ycells = config[0].neighbour_config.y_cells
    var count_zcells = config[0].neighbour_config.z_cells
    var n_cells = count_xcells * count_ycells * count_zcells

    var x_per_super = config[0].neighbour_config.x_cells / config[0].neighbour_config.x_supercells
    var y_per_super = config[0].neighbour_config.y_cells / config[0].neighbour_config.y_supercells
    var z_per_super = config[0].neighbour_config.z_cells / config[0].neighbour_config.z_supercells

    --Compute the bounds for the subcells of this supercell
    var xlo = supercell_id.x * x_per_super
    var xhi = (supercell_id.x + 1) * x_per_super --loops noninclusive so all ok to not -1
    var ylo = supercell_id.y * y_per_super
    var yhi = (supercell_id.y + 1) * y_per_super
    var zlo = supercell_id.z * z_per_super
    var zhi = (supercell_id.z + 1) * z_per_super
for x=xlo, xhi do
for y=ylo, yhi do
for z=zlo, zhi do
    var cell_id : int3d = int3d({x,y,z})
    --Keep track of where to transfer elements and total number of transfers
    var transfer_bounds : int32[27]
    transfer_bounds[0] = 0
    var total_transfers : int32 = 0;
    --Count the number of transfers
    [(function() local __quotes = terralib.newlist()
        for i=1, 26 do
           local queuespace = tradequeues_partitions[i]
           local queue = tradequeues[i]
            __quotes:insert(rquote
                for q in queuespace[cell_id].ispace do
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
    for part in subcell_partition[cell_id].ispace do
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
            local queuespace = tradequeues_partitions[i]
            local queue = tradequeues[i]
            __quotes:insert(rquote
                --Get the bounds for this direction
                var lo = int1d(transfer_bounds[i-1])
                var hi = int1d(transfer_bounds[i])
                for part in subcell_partition[cell_id].ispace do
                    if not(particles[part].neighbour_part_space._valid) then
                        var offset = particles[part].neighbour_part_space._transfer_pos
                        if offset >= lo and offset < hi then
                            var pos : int1d = offset - lo + queuespace[cell_id].bounds.lo;
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
end
end
-- All done
end

--Partitioning code for the tradequeues to correctly create the src/dest tradequeues.
local __demand(__inline) task partition_tradequeue_by_subcells( tradequeue : region(ispace(int1d), part),
                                                                cell_space : ispace(int3d),
                                                                offset : int3d,
                                                                config : region(ispace(int1d), config_type) )
    where reads(config.neighbour_config) do

    var count = tradequeue.bounds.hi + 1
    var count_xcells = config[0].neighbour_config.x_cells
    var count_ycells = config[0].neighbour_config.y_cells
    var count_zcells = config[0].neighbour_config.z_cells
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

local __demand(__inline) task partition_tradequeue_by_supercells( tradequeue : region(ispace(int1d), part),
                                                                  supercell_space : ispace(int3d),
                                                                  offset : int3d,
                                                                  config : region(ispace(int1d), config_type) )
    where reads(config.neighbour_config) do
    var count = tradequeue.bounds.hi + 1
    var count_xcells = config[0].neighbour_config.x_cells
    var count_ycells = config[0].neighbour_config.y_cells
    var count_zcells = config[0].neighbour_config.z_cells
    var n_cells = count_xcells * count_ycells * count_zcells

    var x_per_super = config[0].neighbour_config.x_cells / config[0].neighbour_config.x_supercells
    var y_per_super = config[0].neighbour_config.y_cells / config[0].neighbour_config.y_supercells
    var z_per_super = config[0].neighbour_config.z_cells / config[0].neighbour_config.z_supercells


    --Use legion's coloring option to create this partition - we need a multi domain point coloring since we have multiple domains per color
    var coloring = regentlib.c.legion_multi_domain_point_coloring_create()
    for supercell in supercell_space do
        --The supercell coloring contains all the elements of its subcells, so we loop over the subcells to create it
        var xlo = supercell.x * x_per_super
        var xhi = (supercell.x + 1) * x_per_super --loops noninclusive so all ok to not -1
        var ylo = supercell.y * y_per_super
        var yhi = (supercell.y + 1) * y_per_super
        var zlo = supercell.z * z_per_super
        var zhi = (supercell.z + 1) * z_per_super
    
        --Loop over cells
        for x = xlo, xhi do
            for y = ylo, yhi do
                for z = zlo, zhi do
                    var cell : int3d = int3d({x,y,z})
                    --Ensure we wrap the cell indices (we add # cells to each dimension to handle negative values)
                    var color : int3d = (cell - offset + {count_xcells, count_ycells, count_zcells}) % {count_xcells, count_ycells, count_zcells}
                    --Indices for this color are number of elements per cell (count/n_cells) * (3D to 1D conversion)
                    var oneD_color : int1d = (color.x*count_ycells*count_zcells) + (color.y*count_zcells) + color.z;
                    var rect = rect1d{
                      lo = (count/n_cells)*(oneD_color),
                      hi = (count/n_cells)*(oneD_color+1) - 1
                    }
                    regentlib.c.legion_multi_domain_point_coloring_color_domain(coloring, supercell, rect)
                end
            end
        end
    end
    --Create the partition from the coloring. If offset = 0 this is essentially a controlled equal partition
    var p = partition(disjoint, tradequeue, coloring, supercell_space)
    --Clean up the runtime
    regentlib.c.legion_multi_domain_point_coloring_destroy(coloring)
    return p
end


local __demand(__inline) task partition_sorting_arrays_by_subcells( sort_array : region(ispace(int1d), sorting_ids),
                                                                    particles : region(ispace(int1d), part),
                                                                    cell_space : ispace(int3d),
                                                                    config : region(ispace(int1d), config_type) )
where reads(config), reads(particles.neighbour_part_space.cell_id) do

    
    var count_xcells = config[0].neighbour_config.x_cells
    var count_ycells = config[0].neighbour_config.y_cells
    var count_zcells = config[0].neighbour_config.z_cells
    var num_cells = count_xcells * count_ycells * count_zcells;
    var cell_counts : &int32 = [&int32](regentlib.c.malloc( [terralib.sizeof(int32)] * num_cells))
    var cell_offsets : &int32 = [&int32](regentlib.c.malloc( [terralib.sizeof(int32)] * num_cells))
    var cell_offsets_fixed : &int32 = [&int32](regentlib.c.malloc( [terralib.sizeof(int32)] * num_cells))
    for i = 0, num_cells do
        cell_counts[i] = 0
        cell_offsets[i] = 0
    end
    var copy_space = region(ispace(int1d, 1), part)
    for part in particles do
        --Find the 1D value for the cell this particle belongs to
        var color : int3d = particles[part].neighbour_part_space.cell_id

        var oneD_color : int32 = [int32]((color.x*count_ycells*count_zcells) + (color.y*count_zcells) + color.z);
        cell_counts[oneD_color] = cell_counts[oneD_color] + 1
    end

    --Set offsets
    for i = 1, num_cells do
        cell_offsets[i] = cell_offsets[i-1] + cell_counts[i-1];
    end

    var coloring = regentlib.c.legion_domain_point_coloring_create()
    for color in cell_space do
      --Indices for this color are number of particles (real and padded) in the cell --> (3D to 1D conversion)
      var oneD_color : int32 = int32((color.x*count_ycells*count_zcells) + (color.y*count_zcells) + color.z);
      var rect = rect1d{
        lo = cell_offsets[oneD_color],
        hi = cell_offsets[oneD_color]+cell_counts[oneD_color] - 1
      }
      regentlib.c.legion_domain_point_coloring_color_domain(coloring, color, rect)
    end
    --Create the partition from the coloring. If offset = 0 this is essentially a controlled equal partition
    var p = partition(disjoint, sort_array, coloring, cell_space)
    --Clean up the runtime
    regentlib.c.legion_domain_point_coloring_destroy(coloring)
    return p
end

local __demand(__inline) task partition_sorting_arrays_by_supercells( sort_array : region(ispace(int1d), sorting_ids),
                                                                      particles : region(ispace(int1d), part),
                                                                      supercell_space : ispace(int3d),
                                                                      config : region(ispace(int1d), config_type) )
where reads(config), reads(particles.neighbour_part_space.cell_id) do

    
    var count_xcells = config[0].neighbour_config.x_cells
    var count_ycells = config[0].neighbour_config.y_cells
    var count_zcells = config[0].neighbour_config.z_cells
    var num_cells = count_xcells * count_ycells * count_zcells;
    var cell_counts : &int32 = [&int32](regentlib.c.malloc( [terralib.sizeof(int32)] * num_cells))
    var cell_offsets : &int32 = [&int32](regentlib.c.malloc( [terralib.sizeof(int32)] * num_cells))
    var cell_offsets_fixed : &int32 = [&int32](regentlib.c.malloc( [terralib.sizeof(int32)] * num_cells))
    for i = 0, num_cells do
        cell_counts[i] = 0
        cell_offsets[i] = 0
    end
    var copy_space = region(ispace(int1d, 1), part)
    for part in particles do
        --Find the 1D value for the cell this particle belongs to
        var color : int3d = particles[part].neighbour_part_space.cell_id

        var oneD_color : int32 = [int32]((color.x*count_ycells*count_zcells) + (color.y*count_zcells) + color.z);
        cell_counts[oneD_color] = cell_counts[oneD_color] + 1
    end

    --Set offsets
    for i = 1, num_cells do
        cell_offsets[i] = cell_offsets[i-1] + cell_counts[i-1];
    end

    var x_per_super = config[0].neighbour_config.x_cells / config[0].neighbour_config.x_supercells
    var y_per_super = config[0].neighbour_config.y_cells / config[0].neighbour_config.y_supercells
    var z_per_super = config[0].neighbour_config.z_cells / config[0].neighbour_config.z_supercells


    --Use legion's coloring option to create this partition - we need a multi domain point coloring since we have multiple domains per color
    var coloring = regentlib.c.legion_multi_domain_point_coloring_create()
    for supercell in supercell_space do
        --The supercell coloring contains all the elements of its subcells, so we loop over the subcells to create it
        var xlo = supercell.x * x_per_super
        var xhi = (supercell.x + 1) * x_per_super --loops noninclusive so all ok to not -1
        var ylo = supercell.y * y_per_super
        var yhi = (supercell.y + 1) * y_per_super
        var zlo = supercell.z * z_per_super
        var zhi = (supercell.z + 1) * z_per_super
    
        --Loop over cells
        for x = xlo, xhi do
            for y = ylo, yhi do
                for z = zlo, zhi do
                    var color : int3d = int3d({x,y,z})
                    --Indices for this color are number of elements per cell (count/n_cells) * (3D to 1D conversion)
                    var oneD_color : int32 = int32((color.x*count_ycells*count_zcells) + (color.y*count_zcells) + color.z);
                    var rect = rect1d{
                        lo = cell_offsets[oneD_color],
                        hi = cell_offsets[oneD_color]+cell_counts[oneD_color] - 1
                    }
                    regentlib.c.legion_multi_domain_point_coloring_color_domain(coloring, supercell, rect)
                end
            end
        end
    end
    --Create the partition from the coloring. If offset = 0 this is essentially a controlled equal partition
    var p = partition(disjoint, sort_array, coloring, supercell_space)
    --Clean up the runtime
    regentlib.c.legion_multi_domain_point_coloring_destroy(coloring)
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
        --Avoid Legion issue #1082
        var x = neighbour_init.cell_partition[cell]
        for part in x do
          if x[part].neighbour_part_space._valid then
            regentlib.assert(x[part].neighbour_part_space.cell_id == int3d(cell), "particle found in wrong cell")
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
__demand(__trace)
do
    __demand(__index_launch)
    for slice in [neighbour_init.supercell_partition].colors do
        compute_new_dests( [neighbour_init.supercell_partition][slice], [variables.config]);
    end
--    for cell in [neighbour_init.cell_partition].colors do
--        compute_new_dests( [neighbour_init.cell_partition][cell], [variables.config]);
--    end
    __demand(__index_launch)
    for cell in [neighbour_init.supercell_partition].colors do
        tradequeue_push(cell, [neighbour_init.supercell_partition][cell],
                        [neighbour_init.padded_particle_array], --We pass the full particle array in to reference the partition
                        [neighbour_init.cell_partition], --The subcell partition used to index the computation
                        [variables.config],
                        [generate_range_as_terralist(1, 26):map(function(i) return rexpr
                            [neighbour_init.TradeQueues_bySrc[i]][cell]
                            end
                        end)],
                        [generate_range_as_terralist(1, 26):map(function(i) return rexpr
                            [neighbour_init.TradeQueues[i]]
                            end
                        end)],
                        [generate_range_as_terralist(1, 26):map(function(i) return rexpr
                            [neighbour_init.TradeQueues_bysubSrc[i]]
                            end
                        end)] )
    end
    __demand(__index_launch)
    for cell in [neighbour_init.supercell_partition].colors do
        tradequeue_pull(cell, [neighbour_init.supercell_partition][cell],
                        [neighbour_init.padded_particle_array], --We pass the full particle array in to reference the partition
                        [neighbour_init.cell_partition], --The subcell partition used to index the computation
                        [variables.config],
                        [generate_range_as_terralist(1, 26):map(function(i) return rexpr
                            [neighbour_init.TradeQueues_byDest[i]][cell]
                            end
                        end)],
                        [generate_range_as_terralist(1, 26):map(function(i) return rexpr
                            [neighbour_init.TradeQueues[i]]
                            end
                        end)],
                        [generate_range_as_terralist(1, 26):map(function(i) return rexpr
                            [neighbour_init.TradeQueues_bysubDest[i]]
                            end
                        end)])
    end
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
  var cell_to_supercell_x = [variables.config][0].neighbour_config.x_cells/[variables.config][0].neighbour_config.x_supercells
  var cell_to_supercell_y = [variables.config][0].neighbour_config.y_cells/[variables.config][0].neighbour_config.y_supercells
  var cell_to_supercell_z = [variables.config][0].neighbour_config.z_cells/[variables.config][0].neighbour_config.z_supercells  
  for x=0, [variables.config][0].neighbour_config.x_cells do
    for y=0, [variables.config][0].neighbour_config.y_cells do
      for z=0, [variables.config][0].neighbour_config.z_cells do
         var cell_i3d = int3d({x,y,z})
         var base_index = start_index + (z + y*[variables.config][0].neighbour_config.z_cells + x * [variables.config][0].neighbour_config.y_cells*[variables.config][0].neighbour_config.z_cells)*neighbour_init.padding_per_cell;
        for index = 0, neighbour_init.padding_per_cell do
          --Set the subcell index
          [neighbour_init.padded_particle_array][base_index + index].neighbour_part_space._valid = false;
          [neighbour_init.padded_particle_array][base_index + index].neighbour_part_space.cell_id = cell_i3d;
          --Set the supercell index
          var supercell_id : int3d = int3d({x / cell_to_supercell_x, y / cell_to_supercell_y, z / cell_to_supercell_z});
          [neighbour_init.padded_particle_array][base_index + index].neighbour_part_space.supercell_id = supercell_id;
          [neighbour_init.padded_particle_array][base_index + index].neighbour_part_space.x_cell = int1d(supercell_id.x);

          --Set up the halo indices
          var x_cell_m1 = x - int1d(1)
          if(x_cell_m1 < int1d(0) ) then
              x_cell_m1 = x_cell_m1 + [variables.config][0].neighbour_config.x_cells
          end
          var x_cell_p1 = x + int1d(1)
          if(x_cell_p1 >= int1d([variables.config][0].neighbour_config.x_cells) ) then
              x_cell_p1 = x_cell_p1 - [variables.config][0].neighbour_config.x_cells
          end
          var y_cell_m1 = y - int1d(1)
          if(y_cell_m1 < int1d(0) ) then
              y_cell_m1 = y_cell_m1 + [variables.config][0].neighbour_config.y_cells
          end
          var y_cell_p1 = y + int1d(1)
          if(y_cell_p1 >= int1d([variables.config][0].neighbour_config.y_cells) ) then
              y_cell_p1 = y_cell_p1 - [variables.config][0].neighbour_config.y_cells
          end
          var z_cell_m1 = z - int1d(1)
          if(z_cell_m1 < int1d(0) ) then
              z_cell_m1 = z_cell_m1 + [variables.config][0].neighbour_config.z_cells
          end
          var z_cell_p1 = z + int1d(1)
          if(z_cell_p1 >= int1d([variables.config][0].neighbour_config.z_cells) ) then
              z_cell_p1 = z_cell_p1 - [variables.config][0].neighbour_config.z_cells
          end
          var x_cell = x
          var y_cell = y
          var z_cell = z
         --Save to halos. Note that the "supercell" {-1,-1,-1} doesn't exist and is used to denote that there is no containing halo in this direction
        -- -1, -1, -1
          var temp_supercell : int3d = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z })
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_m1_m1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_m1_m1 = temp_supercell
          end

          -- -1, -1, 0
          temp_supercell = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_m1_0 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_m1_0 = temp_supercell
          end

          -- -1, -1, 1
          temp_supercell = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_m1_p1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_m1_p1 = temp_supercell
          end

          -- -1, 0, -1
          temp_supercell = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_0_m1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_0_m1 = temp_supercell
          end

          -- -1, 0, 0
          temp_supercell = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell / cell_to_supercell_y, z_cell / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_0_0 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_0_0 = temp_supercell
          end

          -- -1, 0, 1
          temp_supercell = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_0_p1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_0_p1 = temp_supercell
          end

          -- -1, 1, -1
          temp_supercell = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_p1_m1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_p1_m1 = temp_supercell
          end

          -- -1, 1, 0
          temp_supercell = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_p1_0 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_p1_0 = temp_supercell
          end

          -- -1, 1, 1
          temp_supercell = int3d({ x_cell_m1 / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_p1_p1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_m1_p1_p1 = temp_supercell
          end

          --0, -1, -1
          temp_supercell = int3d({ x_cell / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_0_m1_m1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_0_m1_m1 = temp_supercell
          end

          -- 0, -1, 0
          temp_supercell = int3d({ x_cell / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_0_m1_0 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_0_m1_0 = temp_supercell
          end

          -- 0, -1, 1
          temp_supercell = int3d({ x_cell / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_0_m1_p1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_0_m1_p1 = temp_supercell
          end

          -- 0, 0, -1
          temp_supercell = int3d({ x_cell / cell_to_supercell_x, y_cell / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_0_0_m1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_0_0_m1 = temp_supercell
          end

          -- 0, 0, 1
          temp_supercell = int3d({ x_cell / cell_to_supercell_x, y_cell / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_0_0_p1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_0_0_p1 = temp_supercell
          end

          -- 0, 1, -1
          temp_supercell = int3d({ x_cell / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_0_p1_m1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_0_p1_m1 = temp_supercell
          end

          -- 0, 1, 0
          temp_supercell = int3d({ x_cell / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_0_p1_0 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_0_p1_0 = temp_supercell
          end

          --0, 1, 1
          temp_supercell = int3d({ x_cell / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_0_p1_p1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_0_p1_p1 = temp_supercell
          end

          --1, -1, -1
          temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_m1_m1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_m1_m1 = temp_supercell
          end

          -- 1, -1, 0
          temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_m1_0 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_m1_0 = temp_supercell
          end

          -- 1, -1, 1
          temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell_m1 / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_m1_p1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_m1_p1 = temp_supercell
          end

          -- 1, 0, -1
          temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_0_m1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_0_m1 = temp_supercell
          end

          -- 1, 0, 0
          temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell / cell_to_supercell_y, z_cell / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_0_0 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_0_0 = temp_supercell
          end

          -- 1, 0, 1
          temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_0_p1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_0_p1 = temp_supercell
          end

          -- 1, 1, -1
          temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell_m1 / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_p1_m1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_p1_m1 = temp_supercell
          end

          -- 1, 1, 0
          temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_p1_0 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_p1_0 = temp_supercell
          end

          -- 1, 1, 1
          temp_supercell = int3d({ x_cell_p1 / cell_to_supercell_x, y_cell_p1 / cell_to_supercell_y, z_cell_p1 / cell_to_supercell_z } )
          if temp_supercell == [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.supercell_id then
              --Not in a halo
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_p1_p1 = int3d({-1, -1, -1})
          else
              [neighbour_init.padded_particle_array][base_index+index].neighbour_part_space.halos_supercell_p1_p1_p1 = temp_supercell
          end

          --All done with halo values 

        end
      end
    end
  end
  --We clear out the originally allocated memory because we don't really want that to exist.
  __delete([variables.particle_array]);

--TODO: Sort particles so that all particles for each cell are contiguous in the region?
--Its possible that we just let the mapper do this instead!
--    var count_xcells = config[0].neighbour_config.x_cells
--    var count_ycells = config[0].neighbour_config.y_cells
--    var count_zcells = config[0].neighbour_config.z_cells
--    var cell_counts : &int32 = [&int32](regentlib.c.malloc( [terralib.sizeof(int32)] * num_cells))
--    var cell_offsets : &int32 = [&int32](regentlib.c.malloc( [terralib.sizeof(int32)] * num_cells))
--    var cell_offsets_fixed : &int32 = [&int32](regentlib.c.malloc( [terralib.sizeof(int32)] * num_cells))
--    for i = 0, num_cells do
--        cell_counts[i] = 0
--        cell_offsets[i] = 0
--    end
--    var copy_space = region(ispace(int1d, 1), part)
--    for part in [neighbour_init.padded_particle_array] do
--        --Find the 1D value for the cell this particle belongs to
--        var color : int3d = [neighbour_init.padded_particle_array][part].neighbour_part_space.cell_id
--
--        var oneD_color : int32 = [int32]((color.x*count_ycells*count_zcells) + (color.y*count_zcells) + color.z);
--        cell_counts[oneD_color] = cell_counts[oneD_color] + 1
--    end
--
--    --Set offsets
--    for i = 1, num_cells do
--        cell_offsets[i] = cell_offsets[i-1] + cell_counts[i-1];
--        cell_offsets_fixed[i] = cell_offsets[i]
--    end
--
--    __delete(copy_space)

--Init Tradequeues
  --The wrapper here is an inline way to create this for all 26 directions inside Regent code, since we're using
  -- symbols to generate this. We could do this with neighbour_init.TradeQueues:map(...) but we may in the future
  --vary tradequeue size based on direction.
  var num_supercells = [variables.config][0].neighbour_config.x_supercells *
                       [variables.config][0].neighbour_config.y_supercells *
                       [variables.config][0].neighbour_config.z_supercells;
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
  --For the source we just divide each of the 26 "direction's" tradequeues into one per source supercell. We use regentlib.coloring
  --code to do this so we know the size/position accurately.
  --For the destination, we do the same, however we use an offset (based upon the direction) to shift each cells tradequeue appropriately.
  var x_cells = [variables.config][0].neighbour_config.x_supercells
  var y_cells = [variables.config][0].neighbour_config.y_supercells
  var z_cells = [variables.config][0].neighbour_config.z_supercells
  var cell_space = ispace(int3d, {x_cells, y_cells, z_cells});
  var x_subcells = [variables.config][0].neighbour_config.x_cells
  var y_subcells = [variables.config][0].neighbour_config.y_cells
  var z_subcells = [variables.config][0].neighbour_config.z_cells
  var cell_space_parameter = ispace(int3d, {x_subcells, y_subcells, z_subcells});
  [(function() local __quotes = terralib.newlist()
    for i = 1, 26 do
        __quotes:insert(rquote
            var [neighbour_init.TradeQueues_bySrc[i]] = partition_tradequeue_by_supercells( [neighbour_init.TradeQueues[i]],
                                                                                       cell_space, int3d({0,0,0}), variables.config );
            var [neighbour_init.TradeQueues_byDest[i]] = partition_tradequeue_by_supercells( [neighbour_init.TradeQueues[i]],
                                                                                        cell_space, [directions[i]], variables.config );
            var [neighbour_init.TradeQueues_bysubSrc[i]] = partition_tradequeue_by_subcells( [neighbour_init.TradeQueues[i]],
                                                                                        cell_space_parameter, int3d({0,0,0}), variables.config );
            var [neighbour_init.TradeQueues_bysubDest[i]] = partition_tradequeue_by_subcells( [neighbour_init.TradeQueues[i]],
                                                                                        cell_space_parameter, [directions[i]], variables.config);
        end)
    end
    return __quotes
   end) ()];

   --Init the supercell partition
    var space_parameter = ispace(int3d, {x_cells, y_cells, z_cells}, {0,0,0})
--    var raw_lp1 = __raw(partition([neighbour_init.padded_particle_array].neighbour_part_space.cell_id, space_parameter));
    var [neighbour_init.supercell_partition] = partition([neighbour_init.padded_particle_array].neighbour_part_space.supercell_id, space_parameter);
    
    var [neighbour_init.cell_partition] = partition([neighbour_init.padded_particle_array].neighbour_part_space.cell_id, cell_space_parameter);

    --Create the halo partition
    var m1_m1_m1_partition  = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_m1_m1_m1, space_parameter)
    var m1_m1_0_partition   = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_m1_m1_0,  space_parameter)
    var m1_m1_p1_partition  = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_m1_m1_p1, space_parameter)
    var m1_0_m1_partition   = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_m1_0_m1,  space_parameter)
    var m1_0_0_partition    = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_m1_0_0,   space_parameter)
    var m1_0_p1_partition   = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_m1_0_p1,  space_parameter)
    var m1_p1_m1_partition  = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_m1_p1_m1, space_parameter)
    var m1_p1_0_partition   = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_m1_p1_0,  space_parameter)
    var m1_p1_p1_partition  = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_m1_p1_p1, space_parameter)
    var _0_m1_m1_partition  = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_0_m1_m1,  space_parameter)
    var _0_m1_0_partition   = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_0_m1_0,   space_parameter)
    var _0_m1_p1_partition  = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_0_m1_p1,  space_parameter)
    var _0_0_m1_partition   = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_0_0_m1,   space_parameter)
    var _0_0_p1_partition   = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_0_0_p1,   space_parameter)
    var _0_p1_m1_partition  = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_0_p1_m1,  space_parameter)
    var _0_p1_0_partition   = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_0_p1_0,   space_parameter)
    var _0_p1_p1_partition  = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_0_p1_p1,  space_parameter)
    var p1_m1_m1_partition  = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_p1_m1_m1, space_parameter)
    var p1_m1_0_partition   = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_p1_m1_0,  space_parameter)
    var p1_m1_p1_partition  = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_p1_m1_p1, space_parameter)
    var p1_0_m1_partition   = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_p1_0_m1,  space_parameter)
    var p1_0_0_partition    = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_p1_0_0,   space_parameter)
    var p1_0_p1_partition   = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_p1_0_p1,  space_parameter)
    var p1_p1_m1_partition  = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_p1_p1_m1, space_parameter)
    var p1_p1_0_partition   = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_p1_p1_0,  space_parameter)
    var p1_p1_p1_partition  = partition([neighbour_init.padded_particle_array].neighbour_part_space.halos_supercell_p1_p1_p1, space_parameter)

    var [neighbour_init.halo_partition] = m1_m1_m1_partition | m1_m1_0_partition  | m1_m1_p1_partition |
                                          m1_0_m1_partition  | m1_0_0_partition   | m1_0_p1_partition  |
                                          m1_p1_m1_partition | m1_p1_0_partition  | m1_p1_p1_partition |
                                          _0_m1_m1_partition | _0_m1_0_partition  | _0_m1_p1_partition |
                                          _0_0_m1_partition  |                      _0_0_p1_partition  |
                                          _0_p1_m1_partition | _0_p1_0_partition  | _0_p1_p1_partition |
                                          p1_m1_m1_partition | p1_m1_0_partition  | p1_m1_p1_partition |
                                          p1_0_m1_partition  | p1_0_0_partition   | p1_0_p1_partition  |
                                          p1_p1_m1_partition | p1_p1_0_partition  | p1_p1_p1_partition
    --Delete all the temporary partitions used to create the halo partition
    __delete(m1_m1_m1_partition)
    __delete(m1_m1_0_partition)
    __delete(m1_m1_p1_partition)
    __delete(m1_0_m1_partition)
    __delete(m1_0_0_partition)
    __delete(m1_0_p1_partition)
    __delete(m1_p1_m1_partition)
    __delete(m1_p1_0_partition)
    __delete(m1_p1_p1_partition)
    __delete(_0_m1_m1_partition)
    __delete(_0_m1_0_partition)
    __delete(_0_m1_p1_partition)
    __delete(_0_0_m1_partition)
    __delete(_0_0_p1_partition)
    __delete(p1_m1_m1_partition)
    __delete(p1_m1_0_partition)
    __delete(p1_m1_p1_partition)
    __delete(p1_0_m1_partition)
    __delete(p1_0_0_partition)
    __delete(p1_0_p1_partition)
    __delete(p1_p1_m1_partition)
    __delete(p1_p1_0_partition)
    __delete(p1_p1_p1_partition)

    --Init the slice partition
    var slice_parameter = ispace(int1d, x_cells, 0);
    var [neighbour_init.x_slices] = partition(complete, [neighbour_init.padded_particle_array].neighbour_part_space.x_cell, slice_parameter);

    --Init the sorting ID arrays
    var [neighbour_init.sorting_id_array] = region(ispace(int1d, tot_parts), sorting_ids);
    --TODO: Divide this into contiguous partitions for cells and non-contiguous partitions for supercells
    var [neighbour_init.sorting_array_cell_partition] = partition_sorting_arrays_by_subcells( [neighbour_init.sorting_id_array],
                                                                                              [neighbour_init.padded_particle_array],
                                                                                              cell_space_parameter,
                                                                                              [variables.config]);

    var [neighbour_init.sorting_array_supercell_partition] = partition_sorting_arrays_by_supercells( [neighbour_init.sorting_id_array],
                                                                                                     [neighbour_init.padded_particle_array],
                                                                                                     cell_space,
                                                                                                     [variables.config]);
end

return initialisation_quote
end

__demand(__inline)
task neighbour_init.check_valid(ns : neighbour_part)

    return ns._valid
    
end

return neighbour_init

-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"
require("defaults")
string_to_field_path = require("src/utils/string_to_fieldpath")

local function recurse_field(field, field_table, type_table, parentstring)
--  for k, v in pairs(field.type) do
    --print(k, v)
--  end
--      print("-------------------------")
--      print(parentstring)
--      print(field.type)
--      for t1, t2 in pairs(field) do
--        print(t1, t2)
--      end
--      print("------")
--      for t1, t2 in pairs(field.type) do
--        print(t1, t2)
--      end
--      print("-------------------------")
  if(field.type.is_fspace_instance ~= nil and field.type.is_fspace_instance) then
    for k, v in pairs(field.type.fields) do
      --print(k, v.field.symbol_name)
      if(parentstring ~= nil) then
        recurse_field(v, field_table, type_table, parentstring.."."..field.field.symbol_name )
      else
        recurse_field(v, field_table,type_table, field.field.symbol_name)
      end
    end
  else
    if(parentstring ~= nil) then
      table.insert(field_table, parentstring.."."..field.field.symbol_name)
      table.insert(type_table, field.field.symbol_type.name)
    else
      table.insert(field_table, field.field.symbol_name)
      table.insert(type_table, field.field.symbol_type.name)
    end
  end
end

local default_value_table = {}
default_value_table["int32"] = rexpr 0 end
default_value_table["uint32"] = rexpr 0 end
default_value_table["int64"] = rexpr 0 end
default_value_table["uint64"] = rexpr 0 end
default_value_table["double"] = rexpr 0.0 end
default_value_table["float"] = rexpr 0.0 end
default_value_table["int1d"] = rexpr int1d(0) end
default_value_table["int2d"] = rexpr int2d({0,0}) end
default_value_table["int3d"] = rexpr int3d({0,0,0}) end

function generate_zero_part_func( )

local field_strings = {}
local string_table = {}
for k, v in pairs(part.fields) do
--  print(v.field.symbol_name)
--  print(v.field.symbol_type)
  recurse_field(v, field_strings, string_table)
end
--print(field_strings)
local mapping_table = terralib.newlist()
for k, _ in pairs(field_strings) do
  if( default_value_table[string_table[k]] == nil) then
    print("No default value set for type: ".. string_table[k]..". Please create an issue to get this added")
  end
  mapping_table:insert({name = string_to_field_path.get_field_path(field_strings[k]), type = string_table[k], default_val = default_value_table[string_table[k]]})
end



local task zero_part_task(particle_region:region(ispace(int1d), part)) where writes(particle_region) do

  for i in particle_region.ispace do
  [mapping_table:map( function(element)
--     print(element.default_val)
     print(element.name)
     return rquote
--     fill(particle_region.[element.name], [element.default_val])
--       fill(particle_region.[element.name], 0)
     particle_region[i].[element.name] = [element.default_val]
     end
  end)];
  end
end

return zero_part_task
end

task zero_core_part(particle_region : region(ispace(int1d), part)) where writes(particle_region.core_part_space) do
fill(particle_region.core_part_space.pos_x, 0.0)
fill(particle_region.core_part_space.pos_y, 0.0)
fill(particle_region.core_part_space.pos_z, 0.0)
fill(particle_region.core_part_space.vel_x, 0.0)
fill(particle_region.core_part_space.vel_y, 0.0)
fill(particle_region.core_part_space.vel_z, 0.0)
fill(particle_region.core_part_space.mass, 0.0)
fill(particle_region.core_part_space.cutoff, 0.0)
fill(particle_region.core_part_space.id, int1d(0))
end

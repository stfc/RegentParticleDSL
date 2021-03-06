-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"
local string_to_field_path = require("src/utils/string_to_fieldpath")
local recursive_fields = require("src/utils/recursive_fields")

local default_value_table = {}
default_value_table["int"] = rexpr 0 end
default_value_table["int32"] = rexpr 0 end
default_value_table["uint32"] = rexpr 0 end
default_value_table["int64"] = rexpr 0 end
default_value_table["int8"] = rexpr int8(0) end
default_value_table["uint64"] = rexpr 0 end
default_value_table["double"] = rexpr 0.0 end
default_value_table["float"] = rexpr 0.0 end
default_value_table["bool"] = rexpr false end
default_value_table["int1d"] = rexpr int1d(0) end
default_value_table["int2d"] = rexpr int2d({0,0}) end
default_value_table["int3d"] = rexpr int3d({0,0,0}) end

function generate_zero_part_quote( particle_array )

local field_strings = {}
local string_table = {}
for k, v in pairs(part.fields) do
--  print(v.field.symbol_name)
--  print(v.field.symbol_type)
  recursive_fields.recurse_field(v, field_strings, string_table)
end
--print(field_strings)
local mapping_table = terralib.newlist()
for k, v in pairs(field_strings) do
  if( default_value_table[string_table[k]] == nil) then
    if type(string_table[k]) == "boolean" then
        print("Uhoh")
    elseif type(string_table[k]) == "table" then
        local field_type = string_table[k]
        if field_type.dim == 1 then
            if default_value_table[field_type.symbol_type_name] ~= nil then
                local t = nil
                if default_value_table[field_type.symbol_type_name .. "[" .. field_type.size[1] .. "]"] == nil then
                    t = terralib.new(field_type.symbol_type[field_type.size[1]])
                    if field_type.symbol_type_name ~= "int1d" and field_type.symbol_type_name ~= "int2d" and field_type.symbol_type_name ~= "int3d" then
                        for i=1, field_type.size[1] do
                            t[i-1] = 0 --default_value_table[field_type.symbol_type_name]
                        end
                    end
                    default_value_table[field_type.symbol_type_name .. "[" .. field_type.size[1] .. "]"] = rexpr t end
                end
                mapping_table:insert({name = string_to_field_path.get_field_path(field_strings[k]), type = string_table[k].symbol_type_name, 
                                        default_val = default_value_table[field_type.symbol_type_name .. "[" .. field_type.size[1] .. "]"]})
            else
                print("No default value set for type: ".. field_type.symbol_type_name ..". Please create an issue to get this added")
            end
        elseif field_type.dim == 2 then
            if default_value_table[field_type.symbol_type_name] ~= nil then
                local t = nil
                if default_value_table[field_type.symbol_type_name .. "[" .. field_type.size[1] .. "][" .. field_type.size[2] .. "]"] == nil then
                    t = terralib.new(field_type.symbol_type[field_type.size[2]][field_type.size[1]])
                    if field_type.symbol_type_name ~= "int1d" and field_type.symbol_type_name ~= "int2d" and field_type.symbol_type_name ~= "int3d" then
                        for j=1, field_type.size[2] do
                            for i=1, field_type.size[1] do
                                t[i-1][j-1] = 0 --default_value_table[field_type.symbol_type_name]
                            end
                        end
                    end
                    default_value_table[field_type.symbol_type_name .. "[" .. field_type.size[1] .. "][" .. field_type.size[2] .. "]"] = rexpr t end
                end
                mapping_table:insert({name = string_to_field_path.get_field_path(field_strings[k]), type = string_table[k].symbol_type_name, 
                    default_val = default_value_table[field_type.symbol_type_name .. "[" .. field_type.size[1] .. "][" .. field_type.size[2] .. "]"]})
            else
                print("No default value set for type: ".. field_type.symbol_type_name ..". Please create an issue to get this added")
            end
        elseif field_type.dim == 3 then
            print("WARNING: Cannot yet automatically initialise arrays of dimension == 3, field ".. field_strings[k].." is not zeroed."..
                     "Make an issue if this is required.")
        else
            print("WARNING: Cannot yet automatically initialise arrays of dimension > 3, field ".. field_strings[k].." is not zeroed."..
                     "Make an issue if this is required.")
        end
    else
        print("No default value set for type: ".. string_table[k]..". Please create an issue to get this added")
    end
  else
    mapping_table:insert({name = string_to_field_path.get_field_path(field_strings[k]), type = string_table[k], default_val = default_value_table[string_table[k]]})
  end
end

return rquote
    [mapping_table:map( function(element)
       return rquote
       fill(particle_array.[element.name], [element.default_val])
       end
    end)];
  end

end

function generate_zero_part_func( )

local field_strings = {}
local string_table = {}
for k, v in pairs(part.fields) do
--  print(v.field.symbol_name)
--  print(v.field.symbol_type)
  recursive_fields.recurse_field(v, field_strings, string_table)
end
--print(field_strings)
local mapping_table = terralib.newlist()
for k, _ in pairs(field_strings) do
  if( default_value_table[string_table[k]] == nil) then
    if type(string_table[k]) == "boolean" then
        print("Uhoh")
    elseif type(string_table[k]) == "table" then
        local field_type = string_table[k]
        if field_type.dim == 1 then
            if default_value_table[field_type.symbol_type_name] ~= nil then
                local t = nil
                if default_value_table[field_type.symbol_type_name .. "[" .. field_type.size[1] .. "]"] == nil then
                    t = terralib.new(field_type.symbol_type[field_type.size[1]])
                    if field_type.symbol_type_name ~= "int1d" and field_type.symbol_type_name ~= "int2d" and field_type.symbol_type_name ~= "int3d" then
                        for i=1, field_type.size[1] do
                            t[i-1] = 0 --default_value_table[field_type.symbol_type_name]
                        end
                    end
                    default_value_table[field_type.symbol_type_name .. "[" .. field_type.size[1] .. "]"] = rexpr t end
                end
                mapping_table:insert({name = string_to_field_path.get_field_path(field_strings[k]), type = string_table[k].symbol_type_name,
                                        default_val = default_value_table[field_type.symbol_type_name .. "[" .. field_type.size[1] .. "]"]})
            else
                print("No default value set for type: ".. field_type.symbol_type_name ..". Please create an issue to get this added")
            end
        elseif field_type.dim == 2 then
            if default_value_table[field_type.symbol_type_name] ~= nil then
                local t = nil
                if default_value_table[field_type.symbol_type_name .. "[" .. field_type.size[1] .. "][" .. field_type.size[2] .. "]"] == nil then
                    t = terralib.new(field_type.symbol_type[field_type.size[2]][field_type.size[1]])
                    if field_type.symbol_type_name ~= "int1d" and field_type.symbol_type_name ~= "int2d" and field_type.symbol_type_name ~= "int3d" then
                        for j=1, field_type.size[2] do
                            for i=1, field_type.size[1] do
                                t[i-1][j-1] = 0 --default_value_table[field_type.symbol_type_name]
                            end
                        end
                    end
                    default_value_table[field_type.symbol_type_name .. "[" .. field_type.size[1] .. "][" .. field_type.size[2] .. "]"] = rexpr t end
                end
                mapping_table:insert({name = string_to_field_path.get_field_path(field_strings[k]), type = string_table[k].symbol_type_name,
                    default_val = default_value_table[field_type.symbol_type_name .. "[" .. field_type.size[1] .. "][" .. field_type.size[2] .. "]"]})
            else
                print("No default value set for type: ".. field_type.symbol_type_name ..". Please create an issue to get this added")
            end    
        elseif field_type.dim == 3 then
            print("WARNING: Cannot yet automatically initialise arrays of dimension == 3, field ".. field_strings[k].." is not zeroed. Make an issue if this is required.")
        else
            print("WARNING: Cannot yet automatically initialise arrays of dimension > 3, field ".. field_strings[k].." is not zeroed. Make an issue if this is required.")
        end
    else
        print("No default value set for type: ".. string_table[k]..". Please create an issue to get this added")
    end
  else
    mapping_table:insert({name = string_to_field_path.get_field_path(field_strings[k]), type = string_table[k], default_val = default_value_table[string_table[k]]})
  end
end



local task zero_part_task(particle_region:region(ispace(int1d), part)) where writes(particle_region) do

--  for i in particle_region.ispace do
  [mapping_table:map( function(element)
     return rquote
     fill(particle_region.[element.name], [element.default_val])
     end
  end)];
--  end
end

return zero_part_task
end

--task zero_core_part(particle_region : region(ispace(int1d), part)) where writes(particle_region.core_part_space) do
--fill(particle_region.core_part_space.pos_x, 0.0)
--fill(particle_region.core_part_space.pos_y, 0.0)
--fill(particle_region.core_part_space.pos_z, 0.0)
--fill(particle_region.core_part_space.vel_x, 0.0)
--fill(particle_region.core_part_space.vel_y, 0.0)
--fill(particle_region.core_part_space.vel_z, 0.0)
--fill(particle_region.core_part_space.mass, 0.0)
--fill(particle_region.core_part_space.cutoff, 0.0)
--fill(particle_region.core_part_space.id, int1d(0))
--end

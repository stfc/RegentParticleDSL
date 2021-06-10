import "regent"
local string_to_field_path = require("src/utils/string_to_fieldpath")
local recursive_fields = require("src/utils/recursive_fields")

local default_value_table = {}
default_value_table["int"] = rexpr 0 end
default_value_table["int8"] = rexpr int8(0) end
default_value_table["int32"] = rexpr 0 end
default_value_table["uint32"] = rexpr 0 end
default_value_table["int64"] = rexpr 0 end
default_value_table["uint64"] = rexpr 0 end
default_value_table["double"] = rexpr 0.0 end
default_value_table["float"] = rexpr 0.0 end
default_value_table["bool"] = rexpr false end
default_value_table["int1d"] = rexpr int1d(0) end
default_value_table["int2d"] = rexpr int2d({0,0}) end
default_value_table["int3d"] = rexpr int3d({0,0,0}) end

local zero_config = {}

function zero_config.generate_zero_config_quote( config ) 
local field_strings = {}
local string_table = {}
for k, v in pairs(config_type.fields) do
  recursive_fields.recurse_field(v, field_strings, string_table)
end

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

local task zero_config_task(config : region(ispace(int1d), config_type)) where writes(config) do
  [mapping_table:map( function(element)
     return rquote
     fill(config.[element.name], [element.default_val])
     end
  end)];
end
    return zero_config_task
end

return zero_config

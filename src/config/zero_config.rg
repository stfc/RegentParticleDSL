import "regent"
local string_to_field_path = require("src/utils/string_to_fieldpath")
local recursive_fields = require("src/utils/recursive_fields")

local default_value_table = {}
default_value_table["int"] = rexpr 0 end
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
--  print(v.field.symbol_name)
--  print(v.field.symbol_type)
  recursive_fields.recurse_field(v, field_strings, string_table)
end

local mapping_table = terralib.newlist()
for k, _ in pairs(field_strings) do
--    print(k)
--    print(field_strings[k])
--    print(string_table[k])
  if( default_value_table[string_table[k]] == nil) then
    if type(string_table[k]) ~= "boolean" then
        print("No default value set for type: ".. string_table[k]..". Please create an issue to get this added")
    else
        print("WARNING: Cannot yet initialise arrays, field ".. field_strings[k].." is not zeroed") 
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

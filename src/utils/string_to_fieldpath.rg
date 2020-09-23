import "regent"

string_to_field_path = {}


local function split_into_table( field_path_string )
  local sep, fields = ".", {}
  local pattern = string.format("([^%s]+)", sep)
  field_path_string:gsub(pattern, function(c) fields[#fields+1] = c end)
  return fields
end

function string_to_field_path.get_field_path( input_string )
  return regentlib.field_path(unpack(split_into_table( input_string ) ) )
end


return string_to_field_path

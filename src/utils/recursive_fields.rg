import "regent"

local recursive_fields = {}

function recursive_fields.recurse_field(field, field_table, type_table, parentstring)
  --If the field is a field space, then we need to recurse until we find a non field space entry
  --This fields name (and any parents) will be prepended in the correct order to the field we find
  if(field.type.is_fspace_instance ~= nil and field.type.is_fspace_instance) then
    for k, v in pairs(field.type.fields) do
      if(parentstring ~= nil) then
        recursive_fields.recurse_field(v, field_table, type_table, parentstring.."."..field.field.symbol_name )
      else
        recursive_fields.recurse_field(v, field_table,type_table, field.field.symbol_name)
      end
    end
  else
     --Our field is not a field space, so we pul lthe name and type and add them to the relevant tables
    if(parentstring ~= nil) then
      table.insert(field_table, parentstring.."."..field.field.symbol_name)
      table.insert(type_table, field.field.symbol_type.name)
    else
      table.insert(field_table, field.field.symbol_name)
      table.insert(type_table, field.field.symbol_type.name)
    end
  end
end

return recursive_fields

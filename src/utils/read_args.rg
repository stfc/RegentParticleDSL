import "regent"

local get_args = {}

function get_args.getarg(key)
  for i=0, #arg-1 do
    if(arg[i] == key) then
      return arg[i+1]
    end
  end
  print("FAILURE: Failed to find input argument")
  os.exit(1)
  return nil
end

function get_args.get_optional_arg(key)
  for i=0, #arg-1 do
    if(arg[i] == key) then
      return arg[i+1]
    end
  end
  return nil
end

return get_args

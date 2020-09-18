import "regent"

--require("defaults")
--require("src/particles/init_part")
--require("src/neighbour_search/cell_pair/neighbour_search")
--require("src/io_modules/HDF5/spaces")

format = require("std/format")
local c = regentlib.c

terralib.includepath = ".;"..terralib.includepath
terralib.includepath = "/home/aidan/hdf5/hdf5-1.12.0/hdf5/include/;"..terralib.includepath
local h5lib = terralib.includec("hdf5.h")
--HDF5 wrapper header as terra can't see some "defines" from hdf5.h.
local wrap = terralib.includec("hdf5_wrapper.h")


--fspace part{
--  test1 : int,
--  test2 : double
--}
--
--fspace hdf5_io_space{
--  name1 : int,
--  name2 : double
--}

simple_hdf5_module = {}

wrap.init_wrapper()
local type_mapping = {}
type_mapping["int"] = wrap.WRAP_H5T_STD_I32LE
type_mapping["int32"] = wrap.WRAP_H5T_STD_I32LE
type_mapping["int64"] = wrap.WRAP_H5T_STD_I64LE
type_mapping["double"] = wrap.WRAP_H5T_NATIVE_DOUBLE
type_mapping["float"] = wrap.WRAP_H5T_NATIVE_FLOAT

local hdf5_io_space = {}

--local variables = {}
--variables.particle_array = regentlib.newsymbol("parts")

--local function create_hdf5_file2(filename, particle_array, hdf5_io_space)
--  local io_fields = terralib.newlist()
--  for k, v in pairs(hdf5_io_space.fields) do
--    io_fields:insert({string=v.field.symbol_name, hdftype=type_mapping[v.field.symbol_type.name]})
--  end
--local create_code = rquote
--  wrap.init_wrapper()
--
--  var file_id = h5lib.H5Fcreate(filename, wrap.WRAP_H5F_ACC_TRUNC, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
--  regentlib.assert(file_id > 0, "Failed to create file")
--
--  var h_dims : h5lib.hsize_t[1]
--  h_dims[0] = particle_array.bounds.hi - particle_array.bounds.lo + 1
--  var dim = h5lib.H5Screate_simple(1, h_dims, [&uint64](0))
--  regentlib.assert(dim > 0, "Failed to create dimensions");
--
----    return rquote format.println("{} {}", [field.string], [field.hdftype]) end
--  [io_fields:map(function(field)
--    return rquote
--       var namething = h5lib.H5Dcreate2(file_id, [field.string], [field.hdftype], dim, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT);
--       regentlib.assert(namething > 0, "Failed to create dataset for [field.string]");
--       h5lib.H5Dclose(namething);
--   end 
--   end)];
--  h5lib.H5Fclose(file_id)
--end
--return create_code
--end

--This local function creates the io field space needed for IO from a given mapper.
local function create_io_fspace(mapper)
 --mapper_fields stores the relationship between the field name (from the mapper) and the 
 --type of the mapped field from the part structure.
 local mapper_fields = terralib.newlist()
 --Loop over the key/value pairs in the mapper (io_field -> part_field)
 for k, v in pairs(mapper) do
   --Search for a matching symbol_name in the particle
   for key, index in pairs(part.fields) do
     if index.field.symbol_name == v then
       --If we find it then we store the io_field name and the type
       mapper_fields:insert({field_name=k, field_type=index.field.symbol_type})
     --We recurse 1 subfield depth so we can see fields inside core_part_space (or other sub field spaces).
     --NB If a field was e.g. field1.field2.field3 we cannot find it.
     elseif string.sub(v, 1, string.len(index.field.symbol_name)) == index.field.symbol_name then
       for key2, index2 in pairs(index.field.symbol_type.fields) do
         if string.sub(v, string.len(index.field.symbol_name)+2) == index2.field.symbol_name then
           mapper_fields:insert({field_name=k, field_type=index2.field.symbol_type})
         end
       end
     end
   end
 end
 local io_type = terralib.types.newstruct("io_type")
 --The idea for how to construct a field space is taken from the regent tests: run_pass/quote_fields.rg
 --We adapted it to allow types to be dependent on the part field space.
 io_type.entries = mapper_fields:map(function(field) return { field.field_name, field.field_type} end)
 return io_type, mapper_fields
end

--This function creates and writes a HDF5 file output, according to the mapper (and corresponding IO_field declaration)
--Arguments: 
--filename : string -- Filename to write to
--mapper : table -- Mapping between IO fields and part type
--particle_array : symbol -- The regent symbol (usually variables.particle_array) representing the particle array.
function simple_hdf5_module.write_output(filename, mapper, particle_array)
  --Create the field space and the io_type_mapping from create_io_fspace
  local hdf5_io_space, io_type_mapping = create_io_fspace(mapper)
  --Mapper fields creates a mapping between the io fspace and part fspace.
  --Due to https://github.com/StanfordLegion/legion/issues/932 we can't just
  --use fields such as core_part_space.pos_x directly, so we check if the part space
  --is a member of core_part_space, and if so store that information and workaround it.
  local mapper_fields = terralib.newlist()
  for k, v in pairs(mapper) do
    local starts_core = (string.sub(v,1,string.len("core_part_space.")) == "core_part_space.")
    local subfield = ""
    if(starts_core) then
      subfield = string.sub(v, string.len("core.part_space.")+1)
    else
      subfield = v
    end
    mapper_fields:insert({io_field=k, part_field=subfield, core_part_field=starts_core})
  end
  --io_fields just stores a list of the fields in the io fspace. We need this for
  --attaching the region correctly
  local io_fields = terralib.newlist()
  for k, v in pairs(mapper) do
    io_fields:insert(k)
  end
  --init_mapping is created from the io_type_mapping, and stores the hdf5 field name (the same as the io space name)
  --and the corresponding hdf5type declaration (taken from the type_mapping)
  local init_mapping = terralib.newlist()
  io_type_mapping:map(function(field) init_mapping:insert({string=field.field_name, hdftype=type_mapping[field.field_type.name] }) end)

  --The copy between the hdf5 region and the particle array needs to be done in a task to work. This is created here
  local task copy_task( hdfreg : region(ispace(int1d), hdf5_io_space), particle_array : region(ispace(int1d), part) ) where
    reads(particle_array), writes(hdfreg) do
    --Loop over the elements of the regions
    for i in hdfreg.ispace do
      --Copy every element from the particle_array to the corresponding field in the io_region.
      --Due to https://github.com/StanfordLegion/legion/issues/932 we check if the field is
      --a member of the core_part_space, and if so use the corresponding branch
      [mapper_fields:map(function(field)
         if(field.core_part_field) then
           return rquote
             hdfreg[i].[field.io_field] = particle_array[i].core_part_space.[field.part_field]
           end
         else
           return rquote
           hdfreg[i].[field.io_field] = particle_array[i].[field.part_field]
           end
         end
      end)];
    end
  
end
  
  local create_code = rquote
    --Create the file
    wrap.init_wrapper()
  
    var file_id = h5lib.H5Fcreate(filename, wrap.WRAP_H5F_ACC_TRUNC, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
    regentlib.assert(file_id > 0, "Failed to create file")
  
    var h_dims : h5lib.hsize_t[1]
    h_dims[0] = particle_array.bounds.hi - particle_array.bounds.lo + 1
    var dim = h5lib.H5Screate_simple(1, h_dims, [&uint64](0))
    regentlib.assert(dim > 0, "Failed to create dimensions");
    --This creates a Dataset for every member of the io mapper. 
    [init_mapping:map(function(field)
      return rquote
         var namething = h5lib.H5Dcreate2(file_id, [field.string], [field.hdftype], dim, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT);
         regentlib.assert(namething > 0, "Failed to create dataset for [field.string]");
         h5lib.H5Dclose(namething);
     end
     end)];
    h5lib.H5Fclose(file_id)
    --Once the file is created, we create a region to attach the file to, and launch the copy_task to write to the file.
    var hdfreg = region(ispace(int1d, particle_array.ispace.bounds.hi - particle_array.ispace.bounds.lo + 1), hdf5_io_space)
    attach(hdf5, hdfreg.[io_fields], filename, regentlib.file_read_write)
    acquire(hdfreg)
    copy_task(hdfreg, particle_array)
    release(hdfreg)
    detach(hdf5, hdfreg)
  end
  return create_code
end

--local test_task = create_hdf5_file2("test.hdf5", variables.particle_array)

--local function create_hdf5_file(filename)
--local create_code = rquote
--  --This must be called at the start of every function involving hdf5 right now
--  wrap.init_wrapper()
--  
--  var file_id = h5lib.H5Fcreate(filename, wrap.WRAP_H5F_ACC_TRUNC, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
--  regentlib.assert(file_id > 0, "Failed to create file")
--
--  var h_dims : h5lib.hsize_t[1]
--  h_dims[0] = 2
--  var dim = h5lib.H5Screate_simple(1, h_dims, [&uint64](0))
--  regentlib.assert(dim > 0, "Failed to create dimensions")
--  
--  var name1thing = h5lib.H5Dcreate2(file_id, "name1", wrap.WRAP_H5T_STD_I32LE, dim, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
--  regentlib.assert(name1thing > 0, "Failed to create dataset for name1")
--  h5lib.H5Dclose(name1thing)
--  
--  var name2thing = h5lib.H5Dcreate2(file_id, "name2", wrap.WRAP_H5T_NATIVE_DOUBLE, dim, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
--  regentlib.assert(name2thing > 0, "Failed to create dataset for name1")
--  h5lib.H5Dclose(name2thing)
--
--  h5lib.H5Fclose(file_id)
--end
--return create_code
--end

--This code is removed until I can work out how to do it/get a response from Regent devs
--local function create_fspace(mapper)
--local  args = terralib.newlist()
--  for k, v in pairs(mapper) do
----    print(k, v)
--     
--    for key, index in pairs(part.fields) do
--      print(key, index)
--      for a, b in pairs(index) do
--        print(a, b)
--      end
--      print("---") 
--      if index.field.symbol_name == v then
----         print(k .." : ".. index.field.symbol_type.name .. ",")
----         args:insert(k .." : ".. index.field.symbol_type.name)
----         print(type(index.field.symbol_type), type(int32))
--         args:insert(regentlib.newsymbol(index.field.symbol_type, k))
----         print(index.field.symbol_type)
----         for key2, index2 in pairs(index.field.symbol_type) do
----           print(key2, index2, type(index2))
----         end
--      end
--    end
--  end
--  print(args)
--  
--local fspace thing{
--
--}
--for key, index in pairs(part) do
--  print(key, index)
--end
--print("---")
--for key, index in pairs(thing) do
--  print(key, index)
--end
--print(args[1])
--return thing
--end

--function simple_hdf5_module.read_file(filename, mapper)
--TODO
--end

--task copy_task(hdfreg : region(ispace(int1d), hdf5_io_space), particle_array : region(ispace(int1d), part)) where
--  reads(particle_array), writes(hdfreg) do
--  for i in particle_array.ispace do
--    hdfreg[i].name1 = particle_array[i].test1
--    hdfreg[i].name2 = particle_array[i].test2
--  end
--end

--mapper = {}
--mapper["name1"] = "test1"
--mapper["name2"] = "test2"


return simple_hdf5_module

--task main()
--  var i0 = ispace(int1d, 2)
--  var [variables.particle_array] = region(i0, part)
--  variables.particle_array[0].test1 = 0
--  variables.particle_array[0].test2 = 1.0
--  variables.particle_array[1].test1 = 1
--  variables.particle_array[1].test2 = 2.0
----  thing([variables.particle_array], "test.hdf5")
--  [write_output("test.hdf5", mapper, variables.particle_array)]
----  [test_task]
----  thing(parts)
--end
----create_fspace(mapper)
--regentlib.start(main)

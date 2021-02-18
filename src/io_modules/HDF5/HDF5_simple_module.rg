-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"


require("src/particles/init_part")

string_to_field_path = require("src/utils/string_to_fieldpath")
format = require("std/format")
local c = regentlib.c

--FIXME: The direct filepath here is bad, should use os.env() instead.
terralib.includepath = ".;"..terralib.includepath
terralib.includepath = os.getenv("HDF5_INCLUDE_PATH")..";"..terralib.includepath
local h5lib = terralib.includec("hdf5.h")
--HDF5 wrapper header as terra can't see some "defines" from hdf5.h.
local wrap = terralib.includec("hdf5_wrapper.h")
local stdlib = terralib.includec("stdlib.h")

simple_hdf5_module = {}

wrap.init_wrapper()
local type_mapping = {}
type_mapping["int"] = wrap.WRAP_H5T_STD_I32LE
type_mapping["int32"] = wrap.WRAP_H5T_STD_I32LE
type_mapping["uint32"] = wrap.WRAP_H5T_STD_U32LE
type_mapping["int64"] = wrap.WRAP_H5T_STD_I64LE
type_mapping["uint64"] = wrap.WRAP_H5T_STD_U64LE
type_mapping["double"] = wrap.WRAP_H5T_NATIVE_DOUBLE
type_mapping["float"] = wrap.WRAP_H5T_NATIVE_FLOAT
type_mapping["int1d"] = wrap.WRAP_H5T_STD_I64LE

local hdf5_io_space = {}

simple_hdf5_module.io_region = regentlib.newsymbol("io_region")

--This local function creates the io field space needed for IO from a given mapper.
local function create_io_fspace(mapper)
 --mapper_fields stores the relationship between the field name (from the mapper) and the 
 --type of the mapped field from the part structure.
 local mapper_fields = terralib.newlist()

 if mapper == nil then
    print("Nil mapper passed into HDF5 IO module, please check a valid mapper is passed into the module")
    os.exit()
 end
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


--This function creates a task that computes how many particles in the file.
--On top of this, it checks that the mapper provided is legal, i.e. that 
--all the fields expected in the file are present.
local function get_particle_count_task(filename, init_mapping)

local task check_file_struc(filename : regentlib.string) 
  --Open the file
  wrap.init_wrapper()

  var file_id = h5lib.H5Fopen(filename, wrap.WRAP_H5F_ACC_RDONLY, wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(file_id >= 0, "Failed to open file");

  [init_mapping:map(function(field)
    return rquote
       var namething = h5lib.H5Dopen2(file_id, [field.string], wrap.WRAP_H5P_DEFAULT)
       if(namething <= 0) then
         var buffer = [rawstring](regentlib.c.malloc(512))
         format.snprintln(buffer, 512, "Failed to find dataset: {}", [field.string])
         regentlib.assert(namething > 0, buffer)
         regentlib.c.free(buffer)
       end
       h5lib.H5Dclose(namething)
   end
   end)];
 
   var firstelem = h5lib.H5Dopen2(file_id, [init_mapping[1].string], wrap.WRAP_H5P_DEFAULT)
   var space = h5lib.H5Dget_space(firstelem)
   var ndims = h5lib.H5Sget_simple_extent_ndims(space)
   regentlib.assert(ndims == 1, "Multi dimensional dataset detected, which we can't handle at this time")
   var dims : h5lib.hsize_t[1]
   regentlib.assert(h5lib.H5Sget_simple_extent_dims(space, dims, [&uint64](0)) > 0, "Failed to read the particle count")
  
   h5lib.H5Dclose(firstelem)
   h5lib.H5Fclose(file_id) 
   format.println("HDF5 read {} parts", dims[0])
   return dims[0]
end
return check_file_struc
end

local zero_part_task = generate_zero_part_func()
--Particle initialisation task - clears out the core_part and neighbour_part values
local task particle_initialisation(particle_region : region(ispace(int1d), part), filename : regentlib.string) where writes(particle_region) do

zero_part_task(particle_region)
end

function simple_hdf5_module.initialise_io_module(particle_array, mapper)
  local hdf5_io_space, io_type_mapping = create_io_fspace(mapper)
  local init_quote = rquote
    var num_parts = [particle_array].ispace.bounds.hi - [particle_array].ispace.bounds.lo + 1
    var io_space = ispace(int1d, num_parts)
    var [simple_hdf5_module.io_region] = region(io_space, hdf5_io_space)
  end
  return init_quote
end

--The initialisation task for the simple hdf5 module
--Creates a particle array and fills it with data from the file.
--Parameters:
--filename: The (relative or absolute) path to the hdf5 input file
--mapper: The mapper table that allows RegentParticleDSL to understand the HDF5 file.
--particle_array: The particle array symbol to be used to store this particle array.
function simple_hdf5_module.read_file(filename, mapper, particle_array)
  --Create the field space and the io_type_mapping from create_io_fspace
  local hdf5_io_space, io_type_mapping = create_io_fspace(mapper)
  --Mapper fields creates a mapping between the io fspace and part fspace.
  --We use the string_to_fieldpath code to allow nested field space accesses
  local mapper_fields = terralib.newlist()
  for k, v in pairs(mapper) do
    mapper_fields:insert({io_field=string_to_field_path.get_field_path(k), part_field=string_to_field_path.get_field_path(v)})
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
  io_type_mapping:map(function(field)
  if type_mapping[field.field_type.name] == nil then
    print("Type "..field.field_type.name.." does not have a known HDF5 conversion type. Please"..
          " open an issue to get this added." )
    os.exit()
  end
  init_mapping:insert({string=field.field_name, hdftype=type_mapping[field.field_type.name] })
  end)

  local get_part_count = get_particle_count_task(filename,  init_mapping)

  --The copy between the hdf5 region and the particle array needs to be done in a task to work. This is created here
  local task copy_task( hdfreg : region(ispace(int1d), hdf5_io_space), particle_array : region(ispace(int1d), part) ) where
    writes(particle_array), reads(hdfreg) do
    --Loop over the elements of the regions
    for i in particle_array.ispace do
      --Copy every element from the particle_array to the corresponding field in the io_region.
      --Due to https://github.com/StanfordLegion/legion/issues/932 we check if the field is
      --a member of the core_part_space, and if so use the corresponding branch
      [mapper_fields:map(function(field)
           return rquote
           particle_array[i].[field.part_field] = hdfreg[i].[field.io_field]
           end
      end)];
    end

end

  local test = rquote
   var part_count = get_part_count(filename)
   format.println("Reading {} particles from the HDF5 file", part_count)
   var particle_space = ispace(int1d, part_count)
   var [particle_array] = region(particle_space, part)

   particle_initialisation([particle_array], filename)

   --Once the other initialisation is done, we read in the particles from the hdf5 file.
   var hdfreg = region(ispace(int1d, part_count), hdf5_io_space)
   attach(hdf5, hdfreg.[io_fields], filename, regentlib.file_read_write)
   acquire(hdfreg)
   copy_task(hdfreg, [particle_array])
   release(hdfreg)
   detach(hdf5, hdfreg)
   c.legion_runtime_issue_execution_fence(__runtime(), __context())
   __delete(hdfreg)
  end
  return test
end


--The initialisation task for the simple hdf5 module
--Creates a particle array and fills it with data from the file.
--Parameters:
--filename: The (relative or absolute) path to the hdf5 input file
--mapper: The mapper table that allows RegentParticleDSL to understand the HDF5 file.
--variables: The variables table that stores the symbols used in the program.
--space_x: The dimension of the x dimension of the box
--space_y: The dimension of the y dimension of the box.
--space_z : Optional. The dimension of the z dimension of the box.
function simple_hdf5_module.initialisation(filename, mapper, variables, space_x, space_y, space_z)
  --Create the field space and the io_type_mapping from create_io_fspace
  local hdf5_io_space, io_type_mapping = create_io_fspace(mapper)
  --Mapper fields creates a mapping between the io fspace and part fspace.
  --We use the string_to_fieldpath code to allow nested field space accesses
  local mapper_fields = terralib.newlist()
  for k, v in pairs(mapper) do
    mapper_fields:insert({io_field=string_to_field_path.get_field_path(k), part_field=string_to_field_path.get_field_path(v)})
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
  io_type_mapping:map(function(field)
  if type_mapping[field.field_type.name] == nil then
    print("Type "..field.field_type.name.." does not have a known HDF5 conversion type. Please"..
          " open an issue to get this added." )
    os.exit()
  end
  init_mapping:insert({string=field.field_name, hdftype=type_mapping[field.field_type.name] })
  end)

  local get_part_count = get_particle_count_task(filename,  init_mapping)

  --Allow 2D input spaces
  if space_z == nil then
    space_z = 0.0
  end
 
  --The copy between the hdf5 region and the particle array needs to be done in a task to work. This is created here
  local task copy_task( hdfreg : region(ispace(int1d), hdf5_io_space), particle_array : region(ispace(int1d), part) ) where
    writes(particle_array), reads(hdfreg) do
    --Loop over the elements of the regions
    for i in particle_array.ispace do
      --Copy every element from the particle_array to the corresponding field in the io_region.
      --Due to https://github.com/StanfordLegion/legion/issues/932 we check if the field is
      --a member of the core_part_space, and if so use the corresponding branch
      [mapper_fields:map(function(field)
           return rquote
           particle_array[i].[field.part_field] = hdfreg[i].[field.io_field]
           end
      end)];
    end

end

  local test = rquote
   regentlib.assert(space_x > 0.0, "X dimension of the space hasn't been set")
   regentlib.assert(space_y > 0.0, "Y dimension of the space hasn't been set")
   var part_count = get_part_count(filename)
   format.println("Reading {} particles from the HDF5 file", part_count)
   var particle_space = ispace(int1d, part_count)
   var [variables.particle_array] = region(particle_space, part)
   var [variables.config] = region(ispace(int1d, 1), config_type)
   fill([variables.config].{space.dim_x, space.dim_y, space.dim_z}, 0.0)
   
   init_space(space_x, space_y, space_z, [variables.config])
   particle_initialisation([variables.particle_array], filename)

   --Once the other initialisation is done, we read in the particles from the hdf5 file.
   var hdfreg = region(ispace(int1d, part_count), hdf5_io_space)
   attach(hdf5, hdfreg.[io_fields], filename, regentlib.file_read_write)
   acquire(hdfreg)
   copy_task(hdfreg, [variables.particle_array])
   release(hdfreg)
   detach(hdf5, hdfreg)
   c.legion_runtime_issue_execution_fence(__runtime(), __context())
  end
  return test
end


function simple_hdf5_module.write_output_manual(filename, mapper, particle_array)
  --Create the field space and the io_type_mapping from create_io_fspace
  local hdf5_io_space, io_type_mapping = create_io_fspace(mapper)
  --Mapper fields creates a mapping between the io fspace and part fspace.
  --We use the string_to_fieldpath code to allow nested field space accesses
  local mapper_fields = terralib.newlist()
  for k, v in pairs(mapper) do
    mapper_fields:insert({io_field=string_to_field_path.get_field_path(k), part_field=string_to_field_path.get_field_path(v)})
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
  io_type_mapping:map(function(field)
  if type_mapping[field.field_type.name] == nil then
    print("Type "..field.field_type.name.." does not have a known HDF5 conversion type. Please"..
          " open an issue to get this added.")
    os.exit()
  end
  init_mapping:insert({string=field.field_name, hdftype=type_mapping[field.field_type.name], regenttype=field.field_type, path=mapper[field.field_name] })
  end)

  local create_code = rquote
    --Create the file
    wrap.init_wrapper()

    var file_id = h5lib.H5Fcreate(filename, wrap.WRAP_H5F_ACC_TRUNC, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
    regentlib.assert(file_id > 0, "Failed to create file")
    var num_parts = 0
    for part in particle_array.ispace do
      var valid = true
      [neighbour_search_validity:map( function(element)
        local field_path = string_to_field_path.get_field_path(element.field)
        return rquote
          valid = valid and ([particle_array][part].[field_path] == [element.result]);
        end
      end)];
      if valid then
        num_parts = num_parts + 1
      end
    end
    var h_dims : h5lib.hsize_t[1]
    h_dims[0] = num_parts
    var dim = h5lib.H5Screate_simple(1, h_dims, [&uint64](0))
    regentlib.assert(dim > 0, "Failed to create dimensions");
    --This creates a Dataset for every member of the io mapper.
    [init_mapping:map(function(field)
      local field_path = string_to_field_path.get_field_path(field.path)
      local stdio = terralib.includec("stdio.h")
      return rquote
         var namething = h5lib.H5Dcreate2(file_id, [field.string], [field.hdftype], dim, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT);
         if(namething <= 0) then
           var buffer = [rawstring](regentlib.c.malloc(512))
           format.snprintln(buffer, 512, "Failed to create dataset: {}", [field.string])
           regentlib.assert(namething > 0, buffer)
           regentlib.c.free(buffer)
         end
         --Once we create the hdf5 dataset store the data for the field
        --FIXME: Array size is not for sure correct...using sizeof(field.regenttype) is not working right now
         var array : &field.regenttype = [&field.regenttype]( c.malloc(8*num_parts))
         var i : int32 = 0
         for part in particle_array do
           var valid : bool = true
           [neighbour_search_validity:map( function(element)
             local sfield_path = string_to_field_path.get_field_path(element.field)
             return rquote
               valid = valid and ([particle_array][part].[sfield_path] == [element.result]);
             end
           end)];
           if valid then
             array[i] = [particle_array][part].[field_path]
             i = i + 1
           end 
         end
         --Now we have the data write it to the file
         var status = h5lib.H5Dwrite(namething, [field.hdftype],  wrap.WRAP_H5S_ALL, wrap.WRAP_H5S_ALL, wrap.WRAP_H5P_DEFAULT, array);
         c.free([array])
         h5lib.H5Dclose(namething);
     end
     end)];
    h5lib.H5Fclose(file_id)
end
return create_code 
end


--This function creates and writes a HDF5 file output, according to the mapper (and corresponding IO_field declaration)
--This function is used from user code.
--Arguments: 
--filename : string -- Filename to write to
--mapper : table -- Mapping between IO fields and part type
--particle_array : symbol -- The regent symbol (usually variables.particle_array) representing the particle array.
function simple_hdf5_module.write_output(filename, mapper, particle_array)
  --Create the field space and the io_type_mapping from create_io_fspace
  local hdf5_io_space, io_type_mapping = create_io_fspace(mapper)
  --Mapper fields creates a mapping between the io fspace and part fspace.
  --We use the string_to_fieldpath code to allow nested field space accesses
  local mapper_fields = terralib.newlist()
  for k, v in pairs(mapper) do
    mapper_fields:insert({io_field=string_to_field_path.get_field_path(k), part_field=string_to_field_path.get_field_path(v)})
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
  io_type_mapping:map(function(field)
  if type_mapping[field.field_type.name] == nil then
    print("Type "..field.field_type.name.." does not have a known HDF5 conversion type. Please"..
          " open an issue to get this added.")
    os.exit()
  end
  init_mapping:insert({string=field.field_name, hdftype=type_mapping[field.field_type.name] })
  end)

  --The copy between the hdf5 region and the particle array needs to be done in a task to work. This is created here
  local task copy_task( hdfreg : region(ispace(int1d), hdf5_io_space), particle_array : region(ispace(int1d), part) ) where
    reads(particle_array), writes(hdfreg) do
    --Loop over the elements of the regions
    var num_written = 0
    for i in particle_array.ispace do
      --Copy every element from the particle_array to the corresponding field in the io_region.
      --Due to https://github.com/StanfordLegion/legion/issues/932 we check if the field is
      --a member of the core_part_space, and if so use the corresponding branch
       var valid = true
      [neighbour_search_validity:map( function(element)
        local field_path = string_to_field_path.get_field_path(element.field)
        return rquote
          valid = valid and (particle_array[i].[field_path] == [element.result]);
        end
       end)];
      if valid then
        [mapper_fields:map(function(field)
             return rquote
             hdfreg[int1d(num_written)].[field.io_field] = particle_array[i].[field.part_field]
             end
        end)];
        num_written = num_written + 1
      end
    end
  
end

  local __demand(__leaf) task count_valid(particles : region(ispace(int1d), part)) : int where reads(particles) do

    var num_parts = 0
    for part in particles.ispace do
      var valid = true
      [neighbour_search_validity:map( function(element)
        local field_path = string_to_field_path.get_field_path(element.field)
        return rquote
          valid = valid and (particles[part].[field_path] == [element.result]);
        end
      end)];
      if valid then
        num_parts = num_parts + 1
      end
    end 
    return num_parts 
  end

  local __demand(__inner) task do_copy(particles : region(ispace(int1d), part), io_reg : region(ispace(int1d), hdf5_io_space), 
                               filename : regentlib.string) where
                               reads(particles, io_reg), writes(particles, io_reg) do
    format.println("{}", filename)
    --Use the premade IO region for IO
    attach(hdf5, io_reg.[io_fields], filename, regentlib.file_read_write)
    acquire(io_reg)
    copy_task(io_reg, particles)
    release(io_reg)
    detach(hdf5, io_reg)
  end
  
  local create_code = rquote
    --Create the file
    wrap.init_wrapper()
    format.println("{}",filename)
    var counter : int = 0
    var file_id = h5lib.H5Fcreate(filename, wrap.WRAP_H5F_ACC_TRUNC, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
    --var file_id = h5lib.H5Fcreate(filename, wrap.WRAP_H5F_ACC_TRUNC, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
    regentlib.assert(file_id > 0, "Failed to create file")
    var num_parts = count_valid(particle_array)
    var h_dims : h5lib.hsize_t[1]
    h_dims[0] = num_parts
    var dim = h5lib.H5Screate_simple(1, h_dims, [&uint64](0))
    regentlib.assert(dim > 0, "Failed to create dimensions");
    --This creates a Dataset for every member of the io mapper. 
    [init_mapping:map(function(field)
      return rquote
         var namething = h5lib.H5Dcreate2(file_id, [field.string], [field.hdftype], dim, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT);
         if(namething <= 0) then
           var buffer = [rawstring](regentlib.c.malloc(512))
           format.snprintln(buffer, 512, "Failed to create dataset: {}", [field.string])
           regentlib.assert(namething > 0, buffer)
           regentlib.c.free(buffer)
         end
         h5lib.H5Dclose(namething);
     end
     end)];
    h5lib.H5Fclose(file_id)
    do_copy(particle_array, [simple_hdf5_module.io_region], filename)
  end
  return create_code
end

return simple_hdf5_module

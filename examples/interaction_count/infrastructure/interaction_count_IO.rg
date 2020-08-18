import "regent"

require("defaults")

local format = require("std/format")
--local h5lib = terralib.includec(os.getenv("HDF_HEADER") or "h5lib.h")
terralib.includepath = ".;"..terralib.includepath
local h5lib = terralib.includec("hdf5.h")
--HDF5 wrapper header as terra can't see some "defines" from hdf5.h.
local wrap = terralib.includec("hdf5_wrapper.h")

stdlib = terralib.includec("stdlib.h")
local c = regentlib.c


terra create_double_array(size : int)
  return stdlib.malloc(sizeof(double) * size)
end

terra create_int_array(size : int)
  return stdlib.malloc(sizeof(int) * size)
end

terra free_double_array(array : &double)
  stdlib.free(array)
end

terra free_int_array(array : &int)
  stdlib.free(array)
end

terra write_field_D(file_id : h5lib.hid_t, group : h5lib.hid_t, fieldname : rawstring, data : &double, size : int)
  --Create the dataset
  var shape : h5lib.hsize_t[1]
  shape[0] = size
  var space = h5lib.H5Screate_simple(1, shape, [&uint64](0))
  var dataset = h5lib.H5Dcreate2( group, fieldname, wrap.WRAP_H5T_IEEE_F64LE, space, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
  
  var status = h5lib.H5Dwrite(dataset,wrap.WRAP_H5T_NATIVE_DOUBLE, wrap.WRAP_H5S_ALL, wrap.WRAP_H5S_ALL, wrap.WRAP_H5P_DEFAULT, data)
  
  h5lib.H5Dclose(dataset)
end

terra write_field_I(file_id: h5lib.hid_t, group : h5lib.hid_t, fieldname : rawstring, data : &int, size : int)
  --Create the dataset
  var shape : h5lib.hsize_t[1]
  shape[0] = size
  var space = h5lib.H5Screate_simple(1, shape, [&uint64](0))
  var dataset = h5lib.H5Dcreate2( group, fieldname, wrap.WRAP_H5T_STD_I32LE, space, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
  
  var status = h5lib.H5Dwrite(dataset,wrap.WRAP_H5T_NATIVE_INT, wrap.WRAP_H5S_ALL, wrap.WRAP_H5S_ALL, wrap.WRAP_H5P_DEFAULT, data)
  
  h5lib.H5Dclose(dataset)

end

task write_hdf5_snapshot(filename : rawstring, particle_list: region(ispace(int1d), part),  config : region(ispace(int1d), config_type)) where reads(particle_list, config) do
  --This must be called at the start of every function involving hdf5 right now
  wrap.init_wrapper()
  
  var file_id = h5lib.H5Fcreate(filename, wrap.WRAP_H5F_ACC_TRUNC, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
  
  --var dims : h5lib.hsize_t[1]
  --format.println("particle_list size {}", particle_list.bounds.hi - particle_list.bounds.lo)
  --dims[0] = particle_list.bounds.hi - particle_list.bounds.lo + 1
  
  --Create file header
  var h_grp = h5lib.H5Gcreate2(file_id, "/Header", wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(h_grp >= 0, "Failed to create file header")
  
  
  --Write the box size
  var h_space = h5lib.H5Screate(wrap.WRAP_H5S_SIMPLE)
  regentlib.assert(h_space >= 0, "failed to create the box size space")
  var dims : h5lib.hsize_t[1]
  dims[0] = 3
  h5lib.H5Sset_extent_simple(h_space, 1, dims, [&uint64](0))
  var h_attr = h5lib.H5Acreate1(h_grp, "BoxSize", wrap.WRAP_H5T_IEEE_F64LE, h_space, wrap.WRAP_H5P_DEFAULT)
  var boxsize : double[3]
  boxsize[0] = config[0].space.dim_x
  boxsize[1] = config[0].space.dim_y
  boxsize[2] = config[0].space.dim_z
  h5lib.H5Awrite(h_attr, wrap.WRAP_H5T_IEEE_F64LE, &boxsize)
  h5lib.H5Sclose(h_space)
  h5lib.H5Aclose(h_attr)
  
  --Close the header
  h5lib.H5Gclose(h_grp)
  
  var part_grp = h5lib.H5Gcreate2(file_id, "/Particles", wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(part_grp >= 0, "Failed to create Particle group")
  --Store particle count as an attribute
  h_space = h5lib.H5Screate(wrap.WRAP_H5S_SIMPLE)
  regentlib.assert(h_space >= 0, "failed to create the particle count space")
  dims[0] = 1
  h5lib.H5Sset_extent_simple(h_space, 1, dims, [&uint64](0))
  h_attr = h5lib.H5Acreate1(part_grp, "ParticleCount", wrap.WRAP_H5T_STD_I64LE, h_space, wrap.WRAP_H5P_DEFAULT)
  var data : int64 = particle_list.bounds.hi - particle_list.bounds.lo + 1
  h5lib.H5Awrite(h_attr, wrap.WRAP_H5T_STD_I64LE, &data)
  h5lib.H5Sclose(h_space)
  h5lib.H5Aclose(h_attr)
  
  --Create the temporary store for the dataset
  var size : int32 = (particle_list.bounds.hi - particle_list.bounds.lo + 1)
  var temp_array : &double = [&double](create_double_array(size)) --[&double](stdlib.malloc(sizeof(double) * size))
  regentlib.assert(temp_array ~= [&double](0), "Failed to allocate the temporary array for TODO:FIELDNAME")
  var counter = 0
  for point in particle_list do
    temp_array[counter] = particle_list[point].core_part_space.pos_x
    counter = counter + 1
  end
  
  write_field_D(file_id, part_grp, "pos_x", temp_array, size)
  
  
  counter = 0
  for point in particle_list do
    temp_array[counter] = particle_list[point].core_part_space.pos_y
    counter = counter + 1
  end
  
  write_field_D(file_id, part_grp, "pos_y", temp_array, size)
  
  
  counter = 0
  for point in particle_list do
    temp_array[counter] = particle_list[point].core_part_space.pos_z
    counter = counter + 1
  end
  
  write_field_D(file_id, part_grp, "pos_z", temp_array, size)
  
  counter = 0
  for point in particle_list do
    temp_array[counter] = particle_list[point].core_part_space.cutoff
    counter = counter + 1
  end
  
  write_field_D(file_id, part_grp, "cutoff", temp_array, size)
  
  free_double_array(temp_array)
  
  var int_array : &int = [&int](create_int_array(size))
  regentlib.assert(int_array ~= [&int](0), "Failed to allocate the temporary array for TODO:FIELDNAME")
  counter = 0
  for point in particle_list do
    int_array[counter] = int(particle_list[point].interactions)
    counter = counter + 1
  end
  
  write_field_I(file_id, part_grp, "interactions", int_array, size)
  
  free_int_array(int_array)
  
  h5lib.H5Gclose(part_grp)
  h5lib.H5Fclose(file_id)

end


task read_particle_count(filename : rawstring) : int64
  --This must be called at the start of every function involving hdf5 right now
  wrap.init_wrapper()
  
  
  var file_id = h5lib.H5Fopen(filename, wrap.WRAP_H5F_ACC_RDONLY, wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(file_id >= 0, "Failed to open file")

  var part_grp = h5lib.H5Gopen2(file_id, "Particles/", wrap.WRAP_H5P_DEFAULT)
  var partcount_attr = h5lib.H5Aopen(part_grp, "ParticleCount", wrap.WRAP_H5P_DEFAULT)
  var part_count : int64
  regentlib.assert( h5lib.H5Aread(partcount_attr, wrap.WRAP_H5T_STD_I64LE, &part_count) >= 0, "Failed to get the ParticleCount")

  h5lib.H5Gclose(part_grp)
  h5lib.H5Fclose(file_id)
  return part_count
end

task read_hdf5_snapshot(filename : rawstring, particle_region : region(ispace(int1d), part) ) where writes(particle_region) do
  --This must be called at the start of every function involving hdf5 right now
  wrap.init_wrapper()
  
  
  var file_id = h5lib.H5Fopen(filename, wrap.WRAP_H5F_ACC_RDONLY, wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(file_id >= 0, "Failed to open file")
  
  var header_grp = h5lib.H5Gopen2(file_id, "Header/", wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(header_grp >= 0, "Failed to open header group")
  var boxsize_attr = h5lib.H5Aopen(header_grp, "BoxSize", wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(boxsize_attr >= 0, "Failed to find the Header/BoxSize")
  
  var boxsize : double[3]
  regentlib.assert( h5lib.H5Aread(boxsize_attr, wrap.WRAP_H5T_IEEE_F64LE, &boxsize) >= 0, "Failed to get the BoxSize attribute")
  format.println("boxsize {} {} {}", boxsize[0], boxsize[1], boxsize[2])
  h5lib.H5Aclose(boxsize_attr)
  h5lib.H5Gclose(header_grp)
  
  var part_grp = h5lib.H5Gopen2(file_id, "Particles/", wrap.WRAP_H5P_DEFAULT)
  var partcount_attr = h5lib.H5Aopen(part_grp, "ParticleCount", wrap.WRAP_H5P_DEFAULT)
  var part_count : int64
  regentlib.assert( h5lib.H5Aread(partcount_attr, wrap.WRAP_H5T_STD_I64LE, &part_count) >= 0, "Failed to get the ParticleCount")
  format.println("Particle Count {}", part_count)
  
  var temp_array : &double = [&double](create_double_array(part_count)) --[&double](stdlib.malloc(sizeof(double) * size))
  regentlib.assert(temp_array ~= [&double](0), "Failed to allocate the temporary array for reading from the file")
  
  var dataset = h5lib.H5Dopen2(part_grp, "pos_x", wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(h5lib.H5Dread(dataset, wrap.WRAP_H5T_NATIVE_DOUBLE, wrap.WRAP_H5S_ALL, wrap.WRAP_H5S_ALL, wrap.WRAP_H5P_DEFAULT, temp_array) >= 0, "Failed to read pos_x array")
  var counter: int = 0
  for point in particle_region do
    particle_region[point].core_part_space.pos_x = temp_array[counter]
    counter = counter + 1
  end
  h5lib.H5Dclose(dataset)
  
  dataset = h5lib.H5Dopen2(part_grp, "pos_y", wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(h5lib.H5Dread(dataset, wrap.WRAP_H5T_NATIVE_DOUBLE, wrap.WRAP_H5S_ALL, wrap.WRAP_H5S_ALL, wrap.WRAP_H5P_DEFAULT, temp_array) >= 0, "Failed to read pos_y array")
  counter= 0
  for point in particle_region do
    particle_region[point].core_part_space.pos_y = temp_array[counter]
    counter = counter + 1
  end
  h5lib.H5Dclose(dataset)
  
  dataset = h5lib.H5Dopen2(part_grp, "pos_z", wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(h5lib.H5Dread(dataset, wrap.WRAP_H5T_NATIVE_DOUBLE, wrap.WRAP_H5S_ALL, wrap.WRAP_H5S_ALL, wrap.WRAP_H5P_DEFAULT, temp_array) >= 0, "Failed to read pos_z array")
  counter= 0
  for point in particle_region do
    particle_region[point].core_part_space.pos_z = temp_array[counter]
    counter = counter + 1
  end
  h5lib.H5Dclose(dataset)
  
  dataset = h5lib.H5Dopen2(part_grp, "cutoff", wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(h5lib.H5Dread(dataset, wrap.WRAP_H5T_NATIVE_DOUBLE, wrap.WRAP_H5S_ALL, wrap.WRAP_H5S_ALL, wrap.WRAP_H5P_DEFAULT, temp_array) >= 0, "Failed to read cutoff array")
  counter= 0
  for point in particle_region do
    particle_region[point].core_part_space.cutoff = temp_array[counter]
    counter = counter + 1
  end
  h5lib.H5Dclose(dataset)
  
  
  h5lib.H5Gclose(part_grp)
  h5lib.H5Fclose(file_id)
end

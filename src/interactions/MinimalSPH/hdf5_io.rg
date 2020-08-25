import "regent"

require("defaults")
require("src/interactions/MinimalSPH/interactions")

local format = require("std/format")
--local h5lib = terralib.includec(os.getenv("HDF_HEADER") or "h5lib.h")
terralib.includepath = ".;"..terralib.includepath
local h5lib = terralib.includec("hdf5.h")
--HDF5 wrapper header as terra can't see some "defines" from hdf5.h.
local wrap = terralib.includec("hdf5_wrapper.h")

stdlib = terralib.includec("stdlib.h")
stdio = terralib.includec("stdio.h")
local c = regentlib.c



function initialisation(filename, particle_array, space)
local initialisation_code = rquote
  var count = read_particle_count(filename)
  format.println("Initialising SPH from {} with {} hydro particles", filename, count)
  var particles_space = ispace(int1d, count)
  var [particle_array] = region(particles_space, part)
  var [space] = region(ispace(int1d, 1), space_config)
  fill([space].{dim_x, dim_y, dim_z}, 0.0)
  read_hdf5_snapshot(filename, count, [particle_array], [space])
end
return initialisation_code
end

terra create_double_array(size : int)
  return stdlib.malloc(sizeof(double) * size)
end

terra create_int_array(size : int)
  return stdlib.malloc(sizeof(int) * size)
end

terra create_float_array(size : int)
  return stdlib.malloc(sizeof(float) * size)
end

terra create_uint64_array(size : int)
  return stdlib.malloc(sizeof(uint64) * size)
end

terra free_double_array(array : &double)
  stdlib.free(array)
end

terra free_int_array(array : &int)
  stdlib.free(array)
end

terra free_float_array( array : &float)
  stdlib.free(array)
end

terra free_uint64_array( array : &uint64)
  stdlib.free(array)
end



task read_particle_count(filename : rawstring) : uint64
  --This must be called at the start of every function involving hdf5 right now
  wrap.init_wrapper()


  var file_id = h5lib.H5Fopen(filename, wrap.WRAP_H5F_ACC_RDONLY, wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(file_id >= 0, "Failed to open file")

  var header_grp = h5lib.H5Gopen2(file_id, "Header/", wrap.WRAP_H5P_DEFAULT)
  var partcount_attr = h5lib.H5Aopen(header_grp, "NumPart_ThisFile", wrap.WRAP_H5P_DEFAULT)
  var part_count : &uint64 = [&uint64](create_uint64_array(6))
  regentlib.assert( h5lib.H5Aread(partcount_attr, wrap.WRAP_H5T_STD_I64LE, part_count) >= 0, "Failed to get the ParticleCount")

  var count = part_count[0]
  free_uint64_array(part_count)
--  format.println("{}", count)
  h5lib.H5Gclose(header_grp)
  h5lib.H5Fclose(file_id)
  return count
end

task read_hdf5_snapshot(filename : rawstring, particle_count : uint64 , particle_region : region(ispace(int1d), part), config_region : region(ispace(int1d), config_type) ) where writes(particle_region, config_region) 
  ,reads(particle_region) do --read access for debugging
  --This must be called at the start of every function involving hdf5 right now
  wrap.init_wrapper()


  var file_id = h5lib.H5Fopen(filename, wrap.WRAP_H5F_ACC_RDONLY, wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(file_id >= 0, "Failed to open file")

--Read the header and boxsize
  var header_grp = h5lib.H5Gopen2(file_id, "Header/", wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(header_grp >= 0, "Failed to open header group")

  var boxsize_attr = h5lib.H5Aopen(header_grp, "BoxSize", wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(boxsize_attr >= 0, "Failed to find the Header/BoxSize")
  var boxsize : double[3]
  regentlib.assert( h5lib.H5Aread(boxsize_attr, wrap.WRAP_H5T_IEEE_F64LE, &boxsize) >= 0, "Failed to get the BoxSize attribute")
  format.println("boxsize {} {} {}", boxsize[0], boxsize[1], boxsize[2])
  config_region[0].space.dim_x = boxsize[0]
  config_region[0].space.dim_y = boxsize[1]
  config_region[0].space.dim_z = boxsize[2]
  h5lib.H5Aclose(boxsize_attr)
  
--Check dimensionality (We can only do 3D tests at the moment
  var dimensionality : int
  var dimension_attr = h5lib.H5Aopen(header_grp, "Dimension", wrap.WRAP_H5P_DEFAULT)
  regentlib.assert( h5lib.H5Aread(dimension_attr, wrap.WRAP_H5T_NATIVE_INT, &dimensionality) >= 0, "Failed to read the Dimension attribute")
  regentlib.assert(dimensionality == 3, "Trying to read a non-3D problem, Aborting")
  h5lib.H5Aclose(dimension_attr)

--Check its only a single file
  var files_per_snapshot : int
  var fps_attr = h5lib.H5Aopen(header_grp, "NumFilesPerSnapshot", wrap.WRAP_H5P_DEFAULT)
  regentlib.assert( h5lib.H5Aread(fps_attr, wrap.WRAP_H5T_NATIVE_INT, &files_per_snapshot) >= 0, "Failed to read the NumFilesPerSnapshot attribute")
  regentlib.assert(files_per_snapshot == 1, "Multiple files per snapshot not currently supported. Aborting")
  h5lib.H5Gclose(header_grp)

--Time to read in the particles
  var parts_grp = h5lib.H5Gopen2(file_id, "PartType0/", wrap.WRAP_H5P_DEFAULT)
  --Lets get the coordinates first
  var coordinates = h5lib.H5Dopen2(parts_grp, "Coordinates", wrap.WRAP_H5P_DEFAULT)
  regentlib.assert( coordinates >= 0, "Failed to open the Coordinates Dataset")
  var coordinate_buffer : &double = [&double](create_double_array(particle_count * 3))
  regentlib.assert(coordinate_buffer ~= [&double](0), "Failed to allocate buffer for the coordinates")
  var shape : h5lib.hsize_t[2]
  var offset : h5lib.hsize_t[2]
  var rank : int
  --Nx3 array
  shape[0] = particle_count
  shape[1] = 3
  offset[0] = 0
  offset[1] = 0
  rank = 2

  --Need to create a memory space so we can read in the shape of the array correctly
  var memspace = h5lib.H5Screate_simple(rank, shape, [&uint64](0))

  --For future safety, we select the (full) hyperslab to read in
  var filespace = h5lib.H5Dget_space(coordinates)
  h5lib.H5Sselect_hyperslab( filespace, wrap.WRAP_H5S_SELECT_SET, offset, [&uint64](0), shape, [&uint64](0))

  --Read the dataspace
  regentlib.assert( h5lib.H5Dread(coordinates, wrap.WRAP_H5T_NATIVE_DOUBLE,  memspace, filespace, wrap.WRAP_H5P_DEFAULT, coordinate_buffer) >= 0,
                   "Failed to read the coordinates")
  var counter = 0
  for part in particle_region.ispace do
    particle_region[part].core_part_space.pos_x = coordinate_buffer[counter*3]
    particle_region[part].core_part_space.pos_y = coordinate_buffer[counter*3+1]
    particle_region[part].core_part_space.pos_z = coordinate_buffer[counter*3+2]
    counter = counter + 1
  end
  --Clean up the coordinate buffer
  h5lib.H5Sclose(filespace)
  h5lib.H5Sclose(memspace)
  h5lib.H5Dclose(coordinates)
  free_double_array(coordinate_buffer)

  --Read the internal energies
  var internal_energy_buffer : &float = [&float](create_float_array(particle_count))
  regentlib.assert(internal_energy_buffer ~= [&float](0), "Failed to allocate the internal energy buffer")
  var internal_energy = h5lib.H5Dopen2(parts_grp, "InternalEnergy", wrap.WRAP_H5P_DEFAULT)
  regentlib.assert( internal_energy >= 0, "Failed to open the InternalEnergy Dataset")
  shape[0] = particle_count
  shape[1] = 1
  offset[0] = 0
  offset[1] = 0
  rank = 1  
  --Need to create a memory space so we can read in the shape of the array correctly
  memspace = h5lib.H5Screate_simple(rank, shape, [&uint64](0))

  --For future safety, we select the (full) hyperslab to read in
  filespace = h5lib.H5Dget_space(internal_energy)
  h5lib.H5Sselect_hyperslab( filespace, wrap.WRAP_H5S_SELECT_SET, offset, [&uint64](0), shape, [&uint64](0))

  regentlib.assert( h5lib.H5Dread(internal_energy, wrap.WRAP_H5T_NATIVE_FLOAT,  memspace, filespace, wrap.WRAP_H5P_DEFAULT, internal_energy_buffer) >= 0,
                   "Failed to read the internal energy")
  counter = 0
  for part in particle_region.ispace do
    particle_region[part].u = internal_energy_buffer[counter]
    counter = counter + 1
  end
  h5lib.H5Sclose(filespace)
  h5lib.H5Sclose(memspace)
  h5lib.H5Dclose(internal_energy)
  free_float_array(internal_energy_buffer)

  var masses_buffer : &float = [&float](create_float_array(particle_count))
  regentlib.assert(masses_buffer ~= [&float](0), "Failed to allocate the masses buffer")
  var masses = h5lib.H5Dopen2(parts_grp, "Masses", wrap.WRAP_H5P_DEFAULT)
  regentlib.assert( masses >= 0, "Failed to open the Masses dataset")
  shape[0] = particle_count
  shape[1] = 1
  offset[0] = 0
  offset[1] = 0
  rank = 1
  --Need to create a memory space so we can read in the shape of the array correctly
  memspace = h5lib.H5Screate_simple(rank, shape, [&uint64](0))

  --For future safety, we select the (full) hyperslab to read in
  filespace = h5lib.H5Dget_space(masses)

  regentlib.assert( h5lib.H5Dread(masses, wrap.WRAP_H5T_NATIVE_FLOAT,  memspace, filespace, wrap.WRAP_H5P_DEFAULT, masses_buffer) >= 0,
                   "Failed to read the masses")
  counter = 0
  for part in particle_region.ispace do
    particle_region[part].core_part_space.mass = masses_buffer[counter]
    counter = counter + 1
  end
  h5lib.H5Sclose(filespace)
  h5lib.H5Sclose(memspace)
  h5lib.H5Dclose(masses)
  free_float_array(masses_buffer)

  var ID_buffer : &uint64 = [&uint64](create_uint64_array(particle_count))
  regentlib.assert(ID_buffer ~= [&uint64](0), "Failed to allocate the ID buffer")
  var IDs = h5lib.H5Dopen2(parts_grp, "ParticleIDs", wrap.WRAP_H5P_DEFAULT)
  regentlib.assert( IDs >= 0, "Failed to open the ParticleIDs dataset")
  shape[0] = particle_count
  shape[1] = 1
  offset[0] = 0
  offset[1] = 0
  rank = 1
  --Need to create a memory space so we can read in the shape of the array correctly
  memspace = h5lib.H5Screate_simple(rank, shape, [&uint64](0))

  --For future safety, we select the (full) hyperslab to read in
  filespace = h5lib.H5Dget_space(IDs)
  regentlib.assert( h5lib.H5Dread(IDs, wrap.WRAP_H5T_NATIVE_ULLONG,  memspace, filespace, wrap.WRAP_H5P_DEFAULT, ID_buffer) >= 0,
                   "Failed to read the particle IDs")
  counter = 0
  for part in particle_region.ispace do
    particle_region[part].core_part_space.id = ID_buffer[counter]
    counter = counter + 1
  end
  h5lib.H5Sclose(filespace)
  h5lib.H5Sclose(memspace)
  h5lib.H5Dclose(IDs)
  free_uint64_array(ID_buffer)

  var smoothing_length_buffer : &float = [&float](create_float_array(particle_count))
  regentlib.assert( smoothing_length_buffer ~= [&float](0), "Failed to allocate the smoothing length buffer")
  var smoothing_lengths = h5lib.H5Dopen2(parts_grp, "SmoothingLength", wrap.WRAP_H5P_DEFAULT)
  regentlib.assert( smoothing_lengths >= 0, "Failed to open the smoothing lengths dataset")
  shape[0] = particle_count
  shape[1] = 1
  offset[0] = 0
  offset[1] = 0
  rank = 1
  --Need to create a memory space so we can read in the shape of the array correctly
  memspace = h5lib.H5Screate_simple(rank, shape, [&uint64](0))

  --For future safety, we select the (full) hyperslab to read in
  filespace = h5lib.H5Dget_space(smoothing_lengths)
  regentlib.assert( h5lib.H5Dread(smoothing_lengths, wrap.WRAP_H5T_NATIVE_FLOAT,  memspace, filespace, wrap.WRAP_H5P_DEFAULT, smoothing_length_buffer) >= 0,
                   "Failed to read the smoothing lengths")

    counter = 0
  for part in particle_region.ispace do
    particle_region[part].h = smoothing_length_buffer[counter]
    --TODO NB. THIS IS A GUESS DO NOT USE PRODUCTION
    particle_region[part].core_part_space.cutoff = kernel_gamma * particle_region[part].h
    counter = counter + 1
  end
  h5lib.H5Sclose(filespace)
  h5lib.H5Sclose(memspace)
  h5lib.H5Dclose(smoothing_lengths)
  free_float_array(smoothing_length_buffer)


  var velocity_buffer : &float = [&float] (create_float_array(particle_count*3))
  regentlib.assert( velocity_buffer ~= [&float](0), "Failed to allocate the velocity buffer")
  var velocities = h5lib.H5Dopen2(parts_grp, "Velocities", wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(velocities >= 0, "Failed to open the velocities dataset")
  --Nx3 array
  shape[0] = particle_count
  shape[1] = 3
  offset[0] = 0
  offset[1] = 0
  rank = 2

  --Need to create a memory space so we can read in the shape of the array correctly
  memspace = h5lib.H5Screate_simple(rank, shape, [&uint64](0))

  --For future safety, we select the (full) hyperslab to read in
  filespace = h5lib.H5Dget_space(velocities)
  regentlib.assert( h5lib.H5Dread(velocities, wrap.WRAP_H5T_NATIVE_FLOAT,  memspace, filespace, wrap.WRAP_H5P_DEFAULT, velocity_buffer) >= 0,
                   "Failed to read the velocities")
  counter = 0
  for part in particle_region.ispace do
    particle_region[part].core_part_space.vel_x = velocity_buffer[counter*3]
    particle_region[part].core_part_space.vel_y = velocity_buffer[counter*3+1]
    particle_region[part].core_part_space.vel_z = velocity_buffer[counter*3+2]
    counter = counter + 1
  end
 
  h5lib.H5Sclose(filespace)
  h5lib.H5Sclose(memspace)
  h5lib.H5Dclose(velocities)
  free_float_array(velocity_buffer)

  --Close the particles group
  h5lib.H5Gclose(parts_grp)
  h5lib.H5Fclose(file_id)
end

terra write_field_D_multidim(file_id : h5lib.hid_t, group : h5lib.hid_t, fieldname : rawstring, data : &double, size : int, dimensions : int)

  var shape : h5lib.hsize_t[2]
  var offsets : h5lib.hsize_t[2]
  var rank : h5lib.hsize_t
  shape[0] = size
  shape[1] = dimensions
  offsets[0] = 0
  offsets[1] = 0
  rank = 2
  var memspace = h5lib.H5Screate(wrap.WRAP_H5S_SIMPLE)
  regentlib.assert( h5lib.H5Sset_extent_simple(memspace, rank, shape, [&uint64](0)) >= 0,"Failed to set the extent")
  var dataset = h5lib.H5Dcreate2(group, fieldname, wrap.WRAP_H5T_IEEE_F64LE, memspace, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(dataset >= 0, "Failed to create the dataset")
  var filespace = h5lib.H5Dget_space(dataset)
  h5lib.H5Sselect_hyperslab(filespace, wrap.WRAP_H5S_SELECT_SET, offsets, [&uint64](0), shape, [&uint64](0) )
  regentlib.assert( h5lib.H5Dwrite( dataset, wrap.WRAP_H5T_NATIVE_DOUBLE, memspace, filespace,  wrap.WRAP_H5P_DEFAULT, data) >= 0, "Failed to write data")
  h5lib.H5Sclose(filespace)
  h5lib.H5Dclose(dataset)
  h5lib.H5Sclose(memspace)
 
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

terra write_field_UL(file_id: h5lib.hid_t, group : h5lib.hid_t, fieldname : rawstring, data : &uint64, size : int)
  --Create the dataset
  var shape : h5lib.hsize_t[1]
  shape[0] = size
  var space = h5lib.H5Screate_simple(1, shape, [&uint64](0))
  var dataset = h5lib.H5Dcreate2( group, fieldname, wrap.WRAP_H5T_NATIVE_ULLONG, space, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)

  var status = h5lib.H5Dwrite(dataset,wrap.WRAP_H5T_NATIVE_ULLONG, wrap.WRAP_H5S_ALL, wrap.WRAP_H5S_ALL, wrap.WRAP_H5P_DEFAULT, data)

  h5lib.H5Dclose(dataset)

end

terra write_field_F_multidim(file_id : h5lib.hid_t, group : h5lib.hid_t, fieldname : rawstring, data : &float, size : int, dimensions : int)
  var shape : h5lib.hsize_t[2]
  var offsets : h5lib.hsize_t[2]
  var rank : h5lib.hsize_t
  shape[0] = size
  shape[1] = dimensions
  offsets[0] = 0
  offsets[1] = 0
  rank = 2
  var memspace = h5lib.H5Screate(wrap.WRAP_H5S_SIMPLE)
  regentlib.assert( h5lib.H5Sset_extent_simple(memspace, rank, shape, [&uint64](0)) >= 0,"Failed to set the extent")
  var dataset = h5lib.H5Dcreate2(group, fieldname, wrap.WRAP_H5T_IEEE_F32LE, memspace, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(dataset >= 0, "Failed to create the dataset")
  var filespace = h5lib.H5Dget_space(dataset)
  h5lib.H5Sselect_hyperslab(filespace, wrap.WRAP_H5S_SELECT_SET, offsets, [&uint64](0), shape, [&uint64](0) )
  regentlib.assert( h5lib.H5Dwrite( dataset, wrap.WRAP_H5T_NATIVE_FLOAT, memspace, filespace,  wrap.WRAP_H5P_DEFAULT, data) >= 0, "Failed to write data")
  h5lib.H5Sclose(filespace)
  h5lib.H5Dclose(dataset)
  h5lib.H5Sclose(memspace)

end


terra write_field_F(file_id: h5lib.hid_t, group : h5lib.hid_t, fieldname : rawstring, data : &float, size : int)
  --Create the dataset
  var shape : h5lib.hsize_t[1]
  shape[0] = size
  var space = h5lib.H5Screate_simple(1, shape, [&uint64](0))
  var dataset = h5lib.H5Dcreate2( group, fieldname, wrap.WRAP_H5T_IEEE_F32LE, space, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)

  var status = h5lib.H5Dwrite(dataset,wrap.WRAP_H5T_NATIVE_FLOAT, wrap.WRAP_H5S_ALL, wrap.WRAP_H5S_ALL, wrap.WRAP_H5P_DEFAULT, data)

  h5lib.H5Dclose(dataset)

end


task write_hdf5_snapshot(filename : rawstring, particle_region : region(ispace(int1d), part), config : region(ispace(int1d), config_type)) where
    reads(particle_region, config) do

  --This must be called at the start of every function involving hdf5 right now
  wrap.init_wrapper()

  var file_id = h5lib.H5Fcreate(filename, wrap.WRAP_H5F_ACC_TRUNC, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(file_id >= 0, "Failed to create file")

  --Create file header
  var h_grp = h5lib.H5Gcreate2(file_id, "/Header", wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(h_grp >= 0, "Failed to create file header")

  --Write the boxsize
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
  h5lib.H5Awrite(h_attr,  wrap.WRAP_H5T_NATIVE_DOUBLE, &boxsize)
  h5lib.H5Sclose(h_space)
  h5lib.H5Aclose(h_attr)
  
  --Write the other attributes that are fixed at the moment
  h_space = h5lib.H5Screate(wrap.WRAP_H5S_SIMPLE)
  regentlib.assert(h_space >= 0, "failed to create the dimension space")
  var dimension = 3
  dims[0] = 1
  h5lib.H5Sset_extent_simple(h_space, 1, dims, [&uint64](0)) 
  h_attr = h5lib.H5Acreate1(h_grp, "Dimension", wrap.WRAP_H5T_NATIVE_INT, h_space, wrap.WRAP_H5P_DEFAULT)
  h5lib.H5Awrite(h_attr, wrap.WRAP_H5T_NATIVE_INT, &dimension)
  h5lib.H5Sclose(h_space)
  h5lib.H5Aclose(h_attr)

  h_space = h5lib.H5Screate(wrap.WRAP_H5S_SIMPLE)
  regentlib.assert(h_space >= 0, "failed to create the entropy ICs space")
  var entropyics = 0
  dims[0] = 1
  h5lib.H5Sset_extent_simple(h_space, 1, dims, [&uint64](0))
  h_attr = h5lib.H5Acreate1(h_grp, "Flag_Entropy_ICs", wrap.WRAP_H5T_STD_I32LE, h_space, wrap.WRAP_H5P_DEFAULT)
  h5lib.H5Awrite(h_attr, wrap.WRAP_H5T_NATIVE_INT, &entropyics)
  h5lib.H5Sclose(h_space)
  h5lib.H5Aclose(h_attr)

  h_space = h5lib.H5Screate(wrap.WRAP_H5S_SIMPLE)
  regentlib.assert(h_space >= 0, "failed to create the MassTable space")
  var masstable : double[6]
  dims[0] = 6
  masstable[0] = 0.0
  masstable[1] = 0.0
  masstable[2] = 0.0
  masstable[3] = 0.0
  masstable[4] = 0.0
  masstable[5] = 0.0
  h5lib.H5Sset_extent_simple(h_space, 1, dims, [&uint64](0))
  h_attr = h5lib.H5Acreate1(h_grp, "MassTable", wrap.WRAP_H5T_IEEE_F64LE, h_space, wrap.WRAP_H5P_DEFAULT)
  h5lib.H5Awrite(h_attr, wrap.WRAP_H5T_NATIVE_DOUBLE, &masstable)
  h5lib.H5Sclose(h_space)
  h5lib.H5Aclose(h_attr)

  h_space = h5lib.H5Screate(wrap.WRAP_H5S_SIMPLE)
  regentlib.assert(h_space >= 0, "failed to create the Files per snapshot space")
  var filepersnap = 1
  dims[0] = 1
  h5lib.H5Sset_extent_simple(h_space, 1, dims, [&uint64](0))
  h_attr = h5lib.H5Acreate1(h_grp, "NumFilesPerSnapshot", wrap.WRAP_H5T_NATIVE_INT, h_space, wrap.WRAP_H5P_DEFAULT)
  h5lib.H5Awrite(h_attr, wrap.WRAP_H5T_NATIVE_INT, &filepersnap)
  h5lib.H5Sclose(h_space)
  h5lib.H5Aclose(h_attr)

  h_space = h5lib.H5Screate(wrap.WRAP_H5S_SIMPLE)
  regentlib.assert(h_space >= 0, "failed to create the NumPartThisFile space")
  var partperfile : int[6]
  dims[0] = 6
  partperfile[0] = particle_region.ispace.bounds.hi - particle_region.ispace.bounds.lo + 1
  partperfile[1] = 0
  partperfile[2] = 0
  partperfile[3] = 0
  partperfile[4] = 0
  partperfile[5] = 0
  h5lib.H5Sset_extent_simple(h_space, 1, dims, [&uint64](0))
  h_attr = h5lib.H5Acreate1(h_grp, "NumPart_ThisFile", wrap.WRAP_H5T_NATIVE_INT, h_space, wrap.WRAP_H5P_DEFAULT)
  h5lib.H5Awrite(h_attr, wrap.WRAP_H5T_NATIVE_INT, &partperfile)
  h5lib.H5Sclose(h_space)
  h5lib.H5Aclose(h_attr)

  h_space = h5lib.H5Screate(wrap.WRAP_H5S_SIMPLE)
  regentlib.assert(h_space >= 0, "failed to create the NumPartTotal space")
  dims[0] = 6
  partperfile[0] = particle_region.ispace.bounds.hi - particle_region.ispace.bounds.lo + 1
  partperfile[1] = 0
  partperfile[2] = 0
  partperfile[3] = 0
  partperfile[4] = 0
  partperfile[5] = 0
  h5lib.H5Sset_extent_simple(h_space, 1, dims, [&uint64](0))
  h_attr = h5lib.H5Acreate1(h_grp, "NumPart_Total", wrap.WRAP_H5T_NATIVE_INT, h_space, wrap.WRAP_H5P_DEFAULT)
  h5lib.H5Awrite(h_attr, wrap.WRAP_H5T_NATIVE_INT, &partperfile)
  h5lib.H5Sclose(h_space)
  h5lib.H5Aclose(h_attr)

  h_space = h5lib.H5Screate(wrap.WRAP_H5S_SIMPLE)
  regentlib.assert(h_space >= 0, "failed to create the NumPartHighWord space")
  dims[0] = 6
  partperfile[0] = 0
  partperfile[1] = 0
  partperfile[2] = 0
  partperfile[3] = 0
  partperfile[4] = 0
  partperfile[5] = 0
  h5lib.H5Sset_extent_simple(h_space, 1, dims, [&uint64](0))
  h_attr = h5lib.H5Acreate1(h_grp, "NumPart_Total_HighWord", wrap.WRAP_H5T_NATIVE_INT, h_space, wrap.WRAP_H5P_DEFAULT)
  h5lib.H5Awrite(h_attr, wrap.WRAP_H5T_NATIVE_INT, &partperfile)
  h5lib.H5Sclose(h_space)
  h5lib.H5Aclose(h_attr)

  h_space = h5lib.H5Screate(wrap.WRAP_H5S_SIMPLE)
  regentlib.assert(h_space >= 0, "failed to create the Time space")
  dims[0] = 1
  --TODO: time will eventually need to be written
  var time = 0.0
  h5lib.H5Sset_extent_simple(h_space, 1, dims, [&uint64](0))
  h_attr = h5lib.H5Acreate1(h_grp, "Time", wrap.WRAP_H5T_NATIVE_DOUBLE, h_space, wrap.WRAP_H5P_DEFAULT)
  h5lib.H5Awrite(h_attr, wrap.WRAP_H5T_NATIVE_DOUBLE, &time)
  h5lib.H5Sclose(h_space)
  h5lib.H5Aclose(h_attr)

  --Close the header
  h5lib.H5Gclose(h_grp)

  --Create the particle group
  var part_grp = h5lib.H5Gcreate2(file_id, "/PartType0", wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT, wrap.WRAP_H5P_DEFAULT)
  regentlib.assert(part_grp >= 0, "Failed to create Particle group")

  --Create the particle group attributes
  h_space = h5lib.H5Screate(wrap.WRAP_H5S_SIMPLE)
  regentlib.assert(h_space >= 0, "failed to create the Number of Fields space")
  var n_fields = 9
  dims[0] = 1
  h5lib.H5Sset_extent_simple(h_space, 1, dims, [&uint64](0))
  h_attr = h5lib.H5Acreate1(part_grp, "NumberOfFields", wrap.WRAP_H5T_NATIVE_INT, h_space, wrap.WRAP_H5P_DEFAULT)
  h5lib.H5Awrite(h_attr, wrap.WRAP_H5T_NATIVE_INT, &n_fields)
  h5lib.H5Sclose(h_space)
  h5lib.H5Aclose(h_attr)

  --Create the particle group attributes
  h_space = h5lib.H5Screate(wrap.WRAP_H5S_SIMPLE)
  regentlib.assert(h_space >= 0, "failed to create the Number of Particles space")
  var n_particles = particle_region.ispace.bounds.hi - particle_region.ispace.bounds.lo + 1 
  dims[0] = 1
  h5lib.H5Sset_extent_simple(h_space, 1, dims, [&uint64](0))
  h_attr = h5lib.H5Acreate1(part_grp, "NumberOfParticles", wrap.WRAP_H5T_NATIVE_INT, h_space, wrap.WRAP_H5P_DEFAULT)
  h5lib.H5Awrite(h_attr, wrap.WRAP_H5T_NATIVE_INT, &n_particles)
  h5lib.H5Sclose(h_space)
  h5lib.H5Aclose(h_attr)

  var part_count = particle_region.ispace.bounds.hi - particle_region.ispace.bounds.lo + 1
  var counter = 0

  --Write the coordinate array
  var coordinate_array : &double = [&double](create_double_array( part_count * 3 ))
  for part in particle_region.ispace do
    coordinate_array[counter*3] = particle_region[part].core_part_space.pos_x
    coordinate_array[counter*3 + 1] = particle_region[part].core_part_space.pos_y
    coordinate_array[counter*3 + 2] = particle_region[part].core_part_space.pos_z
    counter = counter + 1
  end
  write_field_D_multidim(file_id , part_grp , "Coordinates", coordinate_array, part_count, 3)
  free_double_array(coordinate_array)

  var density_array : &float = [&float](create_float_array( part_count ))
  counter = 0
  for part in particle_region.ispace do
    density_array[counter] = particle_region[part].rho
    counter = counter + 1
  end
  write_field_F(file_id, part_grp, "Densities", density_array, part_count)
  free_float_array(density_array)

  var entropy_array : &float = [&float](create_float_array( part_count ))
  counter = 0
  for part in particle_region.ispace do
    entropy_array[counter] = gas_entropy_from_internal_energy(particle_region[part].rho, particle_region[part].u)
    counter = counter + 1
  end
  write_field_F(file_id, part_grp, "Entropies", entropy_array, part_count)
  free_float_array(entropy_array)

  var internal_energy_array : &float = [&float](create_float_array( part_count ))
  counter = 0
  for part in particle_region.ispace do
    internal_energy_array[counter] = particle_region[part].u
    counter = counter + 1
  end
  write_field_F(file_id, part_grp, "InternalEnergies", internal_energy_array, part_count)
  free_float_array(internal_energy_array)

  var masses_array : &float = [&float](create_float_array( part_count ))
  counter = 0
  for part in particle_region.ispace do
    masses_array[counter] = particle_region[part].core_part_space.mass
    counter = counter + 1
  end
  write_field_F(file_id, part_grp, "Masses", masses_array, part_count)
  free_float_array(masses_array)

  var id_array : &uint64 = [&uint64](create_uint64_array( part_count ))
  counter = 0
  for part in particle_region.ispace do
    id_array[counter] = particle_region[part].core_part_space.id
    counter = counter + 1
  end
  write_field_UL(file_id, part_grp, "ParticleIDs", id_array, part_count)
  free_uint64_array(id_array)

  var pressure_array : &float = [&float](create_float_array( part_count ))
  counter = 0
  for part in particle_region.ispace do
    pressure_array[counter] = particle_region[part].pressure
    counter = counter + 1
  end
  write_field_F(file_id, part_grp, "Pressures", pressure_array, part_count)
  free_float_array(pressure_array)

  var smoothing_array : &float = [&float](create_float_array( part_count ))
  counter = 0
  for part in particle_region.ispace do
    smoothing_array[counter] = particle_region[part].h
    counter = counter + 1
  end
  write_field_F(file_id, part_grp, "SmoothingLengths", smoothing_array, part_count)
  free_float_array(smoothing_array)

  var velocity_array : &float = [&float](create_float_array( 3 * part_count ))
  counter = 0
  for part in particle_region.ispace do
    smoothing_array[counter * 3] = particle_region[part].core_part_space.vel_x
    smoothing_array[counter * 3 + 1] = particle_region[part].core_part_space.vel_y
    smoothing_array[counter * 3 + 2] = particle_region[part].core_part_space.vel_z
    counter = counter + 1
  end
  write_field_F_multidim(file_id , part_grp , "Velocities", velocity_array, part_count, 3)
  free_float_array(velocity_array)


  --Close the particle group
  h5lib.H5Gclose(part_grp)

  --Close the file
  h5lib.H5Fclose(file_id)
end




--task main()
--  var count = read_particle_count("/home/aidan/swiftsim/examples/HydroTests/SodShock_3D/sodShock.hdf5")
--  format.println("{}", count)
--  var particles_space = ispace(int1d, count)
--  var particle_region = region(particles_space, part)
--  var space = region(ispace(int1d, 1), space_config)
--  fill(space.{dim_x, dim_y, dim_z}, 0.0)
--  read_hdf5_snapshot("/home/aidan/swiftsim/examples/HydroTests/SodShock_3D/sodShock.hdf5", count, particle_region, space)
--  write_hdf5_snapshot("/home/aidan/test.hdf5", particle_region, space)
--end

--regentlib.start(main)

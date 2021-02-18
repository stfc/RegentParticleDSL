-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"


fspace config_type{
  space : space_config_type,
  neighbour_config : neighbour_config_type,
  timing_config : timing_config_type
  
}

task init_space( dim_x : double, dim_y : double, dim_z : double, config : region(ispace(int1d), config_type)) where
  writes(config) do
  config[0].space.dim_x = dim_x
  config[0].space.dim_y = dim_y
  config[0].space.dim_z = dim_z
  config[0].space.timestep = 0.0
end

task zero_space( config : region(ispace(int1d), config_type) ) where writes(config) do
  init_space( 0.0, 0.0, 0.0, config)
end

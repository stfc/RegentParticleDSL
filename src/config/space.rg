import "regent"

fspace space_config_type{
  dim_x : double, 
  dim_y : double,
  dim_z : double,
  timestep : double
}



task init_space( dim_x : double, dim_y : double, dim_z : double, space : region(ispace(int1d), space_config_type)) where
  writes(space) do
  space[0].dim_x = dim_x 
  space[0].dim_y = dim_y 
  space[0].dim_z = dim_z 
  space[0].timestep = 0.0
end

task zero_space( space : region(ispace(int1d), space_config_type) ) where writes(space) do
  init_space( 0.0, 0.0, 0.0, space)
end

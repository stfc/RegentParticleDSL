import "regent"

fspace space_config{
  dim_x : double, 
  dim_y : double,
  dim_z : double
}


task init_space( dim_x : double, dim_y : double, dim_z : double, space : region(ispace(int1d), space_config)) where
  writes(space) do
  space[0].dim_x = dim_x 
  space[0].dim_y = dim_y 
  space[0].dim_z = dim_z 

end

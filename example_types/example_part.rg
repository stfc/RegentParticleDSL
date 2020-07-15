import "regent"

--This shows an example construction of a custom particle type, used in a computation
--TODO: For now this is hard coded as relative paths are more complex
--This requires TERRA_PATH to be set to the base directory
require("src/particles/core_part")
require("src/neighbour_search/cell_pair/cell_pair")

--This structure should always be called "part"
fspace part{
--We are required to include the neighbour_part and core_part types in our particle
--declaration
  neighbour_part_space : neighbour_part,
  core_part_space : core_part,
--Any extra defitions we want for our computation
  extra_variable_1 : double,
  extra_variable_2 : float,
  extra_variable_3 : int32

}

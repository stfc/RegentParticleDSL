import "regent"

require("src/particles/core_part")

--Settings table
DSL_settings = {}
DSL_settings.DIMENSIONALITY = 3
DSL_settings.PERIODICITY = true
DSL_settings.TIMING = false

--Set empty neighbour init global variable
neighbour_init = {}

--Initialise the variables list
variables = {}

variables.config = regentlib.newsymbol("config")
variables.particle_array = regentlib.newsymbol("particle_region")
variables.io_array = regentlib.newsymbol("io_array")

function set_dimensionality(dimensions)
  if(dimensions ~= 2 and dimensions ~= 3 ) then
    print("Dimensionality other than 2D and 3D not currently supported")
    os.exit(1)
  end
  DSL_settings.DIMENSIONALITY = dimensions
end

function set_periodicity(periodicity)
  if(periodicity ~= true and periodicity ~= false) then
    print("Periodicity input must be a boolean")
    os.exit(1)
  end
  DSL_settings.PERIODICITY = periodicity
end

function enable_timing()
  DSL_settings.TIMING = true
end

function setup_part()
--Periodic imports
  if DSL_settings.PERIODICITY then
    if DSL_settings.DIMENSIONALITY == 3 then
      require("src/neighbour_search/cell_pair_tradequeues/import_cell_pair")
    else
      require("src/neighbour_search/2d_cell_pair_tradequeues/import_cell_pair")
    end
  else
    if DSL_settings.DIMENSIONALITY == 3 then
      require("src/neighbour_search/cell_pair_tradequeues_nonperiod/import_3d_nonperiod")
    else
      require("src/neighbour_search/cell_pair_tradequeues_nonperiod/import_2d_nonperiod")
    end
  end

end


function setup_dsl()
require("src/config/space")
require("src/config/timing")
require("src/config/default_config")
--Periodic imports
if DSL_settings.PERIODICITY then
  if DSL_settings.DIMENSIONALITY == 3 then
    require("src/neighbour_search/cell_pair_tradequeues/neighbour_search")
    neighbour_init = require("src/neighbour_search/cell_pair_tradequeues/neighbour_init")
  else
    require("src/neighbour_search/2d_cell_pair_tradequeues/neighbour_search")
    neighbour_init = require("src/neighbour_search/2d_cell_pair_tradequeues/neighbour_init")
  end
else
  require("src/neighbour_search/cell_pair_tradequeues_nonperiod/neighbour_search")
  neighbour_init = require("src/neighbour_search/cell_pair_tradequeues_nonperiod/neighbour_init")
end




  require("src/utils/invoke_framework")

end

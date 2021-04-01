import "regent"

require("src/particles/core_part")

--Settings table
DSL_settings = {}
DSL_settings.DIMENSIONALITY = 3
DSL_settings.PERIODICITY = true
DSL_settings.TIMING = false
DSL_settings.part_setup = false
DSL_settings.dsl_setup = false

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


--Mark this as completed
  DSL_settings.part_setup = true
end


function setup_dsl()
require("src/config/space")
require("src/config/timing")
require("src/config/default_config")

--Check part setup was done
if DSL_settings.part_setup == false then
    print("Particle setup must be completed before the setup_dsl call")
    os.exit(1)
end


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
  --Mark DSL as setup
   DSL_settings.dsl_setup = true

end

--Lets run the DSL immediately. This is currently just a wrapper for regentlib.start, but lets the DSL
--make mapper choices etc. depending on chosen options/kernels in the future.
function run_DSL( main_function )
  if DSL_settings.dsl_setup == false then
    print("DSL setup must be completed before the run_dsl call")
    os.exit(1)
  end
  regentlib.start(main_function)
end


local terra set_mappers()

end

--Single function for compiling the DSL
function compile_DSL( main_function, executable_name )

  if DSL_settings.dsl_setup == false then
    print("DSL setup must be completed before the compile_dsl call")
    os.exit(1)
  end
  if executable_name == nil then
    executable_name = "a.out"
    print("WARNING: Executable name was not given, naming it a.out")
  end
  --Build function
  local root_dir = "./"
  local out_dir = (os.getenv('OBJNAME') and os.getenv('OBJNAME'):match('.*/')) or root_dir
  local link_flags = terralib.newlist({"-L" .. out_dir, "-lm", "-lhdf5"})
  local exe = os.getenv('OBJNAME') or executable_name
  regentlib.saveobj(main_task, exe, "executable", set_mappers, link_flags)
end

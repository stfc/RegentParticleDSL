import "regent"

require("src/particles/core_part")

--Settings table
DSL_settings = {}
DSL_settings.DIMENSIONALITY = 3
DSL_settings.PERIODICITY = true
DSL_settings.TIMING = false
DSL_settings.part_setup = false
DSL_settings.dsl_setup = false
DSL_settings.ALLTOALL = false
DSL_settings.HIGHPERFORMANCE = true

DSL_settings.mapper_path = nil

--Set empty neighbour init global variable
neighbour_init = {}

--Initialise the variables list
variables = {}

variables.config = regentlib.newsymbol("config")
variables.particle_array = regentlib.newsymbol("particle_region")
variables.io_array = regentlib.newsymbol("io_array")

function disable_high_performance()
    DSL_settings.HIGHPERFORMANCE = false
end

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

--NB: Only supports 3D periodic
function enable_all_to_all()
  DSL_settings.ALLTOALL = true
end

function setup_part()
--Periodic imports
  if DSL_settings.PERIODICITY then
    if DSL_settings.DIMENSIONALITY == 3 then
        if DSL_settings.ALLTOALL then
            require("src/neighbour_search/all_to_all/import_alltoall")
            DSL_settings.mapper_path = "src/mappers/tradequeue_mapper.cc"
        else
            if DSL_settings.HIGHPERFORMANCE then
                require("src/neighbour_search/HP_cell_pair_tradequeues/import_cell_pair")
                DSL_settings.mapper_path = "src/mappers/high_perform_mapper.cc"
            else
                require("src/neighbour_search/cell_pair_tradequeues/import_cell_pair")
                DSL_settings.mapper_path = "src/mappers/tradequeue_mapper.cc"
            end
        end
    else
      require("src/neighbour_search/2d_cell_pair_tradequeues/import_cell_pair")
      DSL_settings.mapper_path = "src/mappers/tradequeue_mapper.cc"
    end
  else
    if DSL_settings.DIMENSIONALITY == 3 then
      require("src/neighbour_search/cell_pair_tradequeues_nonperiod/import_3d_nonperiod")
      DSL_settings.mapper_path = "src/mappers/tradequeue_mapper.cc"
    else
      require("src/neighbour_search/cell_pair_tradequeues_nonperiod/import_2d_nonperiod")
      DSL_settings.mapper_path = "src/mappers/tradequeue_mapper.cc"
    end
  end


--Mark this as completed
  DSL_settings.part_setup = true
end

dl_meso_config_mod = nil
dl_meso_field_mod = nil
dl_meso_start_mod = nil
dl_meso_statistics_mod = nil
dl_meso_timing_mod = nil
dl_meso_init = nil
dl_meso_read_mod = nil
dl_meso_write_mod = nil

function import_dl_meso( custom_config_path, custom_part_path)
    setup_part()
    if custom_part_path ~= nil then
        require(custom_part_path)
    else
        require("src/io_modules/dl_meso/dl_meso_particle")
    end
    require("src/config/space")
    require("src/config/timing")
    if custom_config_path ~= nil then
        require(custom_config_path)
    else
        require("src/io_modules/dl_meso/dl_meso_config")
    end

    if DSL_settings.PERIODICITY then
      if DSL_settings.DIMENSIONALITY == 3 then
        if DSL_settings.ALLTOALL then
            require("src/neighbour_search/all_to_all/neighbour_search")
            neighbour_init = require("src/neighbour_search/all_to_all/neighbour_init")
        else
            if DSL_settings.HIGHPERFORMANCE then
                require("src/neighbour_search/HP_cell_pair_tradequeues/neighbour_search")
                neighbour_init = require("src/neighbour_search/HP_cell_pair_tradequeues/neighbour_init")
            else
                require("src/neighbour_search/cell_pair_tradequeues/neighbour_search")
                neighbour_init = require("src/neighbour_search/cell_pair_tradequeues/neighbour_init")
            end
        end
      else
        require("src/neighbour_search/2d_cell_pair_tradequeues/neighbour_search")
        neighbour_init = require("src/neighbour_search/2d_cell_pair_tradequeues/neighbour_init")
      end
    else
      require("src/neighbour_search/cell_pair_tradequeues_nonperiod/neighbour_search")
      neighbour_init = require("src/neighbour_search/cell_pair_tradequeues_nonperiod/neighbour_init")
    end
    require("src/utils/invoke_framework")
    DSL_settings.dsl_setup = true
    dl_meso_read_mod = require("src/io_modules/dl_meso/read_dl_meso")
    dl_meso_write_mod = require("src/io_modules/dl_meso/write_dl_meso")
    dl_meso_config_mod = require("src/io_modules/dl_meso/dl_meso_config_mod")
    dl_meso_field_mod = require("src/io_modules/dl_meso/dl_meso_field_mod")
    dl_meso_start_mod = require("src/io_modules/dl_meso/dl_meso_start_mod")
    dl_meso_statistics_mod = require("src/io_modules/dl_meso/dl_meso_statistics_mod")
    dl_meso_timing_mod = require("src/io_modules/dl_meso/dl_meso_timing")
    dl_meso_init = require("src/io_modules/dl_meso/dl_meso_initialise")
    print("Imported DL_MESO modules correctly")
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
        if DSL_settings.ALLTOALL then
            require("src/neighbour_search/all_to_all/neighbour_search")
            neighbour_init = require("src/neighbour_search/all_to_all/neighbour_init")
        else
            if DSL_settings.HIGHPERFORMANCE then
                require("src/neighbour_search/HP_cell_pair_tradequeues/neighbour_search")
                neighbour_init = require("src/neighbour_search/HP_cell_pair_tradequeues/neighbour_init")
            else
                require("src/neighbour_search/cell_pair_tradequeues/neighbour_search")
                neighbour_init = require("src/neighbour_search/cell_pair_tradequeues/neighbour_init")
            end
        end
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

function test_mapper(mapper_path)
  DSL_settings.mapper_path = mapper_path
  compile_mapper_run()
end

function compile_mapper_run()
    if DSL_settings.mapper_path == nil then
        return
    end
    --Compile and link the mapper.
    --Find mapper file name
    local fileindices = string.find(DSL_settings.mapper_path, "[A-Za-z0-9_%-]*.cc")
    local filename = string.sub(DSL_settings.mapper_path, fileindices)
    local includepath = string.sub(DSL_settings.mapper_path, 0, fileindices-1)

    --Setup includepath for compiler
    local root_dir = "./"
    local include_path = ""
    local include_dirs = terralib.newlist()
    include_dirs:insert("-I")
    include_dirs:insert(root_dir)
    for path in string.gmatch(os.getenv("INCLUDE_PATH"), "[^;]+") do
        include_path = include_path .. " -I " .. path
        include_dirs:insert("-I")
        include_dirs:insert(path)
    end
    include_dirs:insert("-I")
    include_dirs:insert(includepath)
    include_path = include_path .. " -I " .. includepath

    local mapper_cc = DSL_settings.mapper_path
    local mapper_so = "./libmapper" .. ".so"
    local cxx = os.getenv('CXX') or 'c++'
    local cxx_flags = os.getenv('CXXFLAGS') or ''
    cxx_flags = cxx_flags .. " -g -O2 -Wall -Werror"
    --Not yet working on Mac. circuit_mapper has some stuff to detect Darwin
    cxx_flags = cxx_flags .. " -shared -fPIC"
    local cmd = (cxx .. " " .. cxx_flags .. " " .. include_path .. " " ..
                mapper_cc .. " -o " .. mapper_so)
    if os.execute(cmd) ~= 0 then
        print("Error: failed to compile " .. mapper_cc)
        assert(false)
    end

    regentlib.linklibrary(mapper_so)
    local cmapper = terralib.includec(includepath .. string.gsub(filename, "%.cc", ".h"), include_dirs)
    return cmapper, mapper_so
    
    
end

--Lets run the DSL immediately. This is currently just a wrapper for regentlib.start, but lets the DSL
--make mapper choices etc. depending on chosen options/kernels in the future.
function run_DSL( main_function )
  if DSL_settings.dsl_setup == false then
    print("DSL setup must be completed before the run_dsl call")
    os.exit(1)
  end
  print(DSL_settings.mapper_path)
  if DSL_settings.mapper_path ~= nil then
    local cmapper, _ = compile_mapper_run()
    regentlib.start(main_function, cmapper.register_mappers)
    --regentlib.start(main_function, mapper function TODO
  else
    regentlib.start(main_function)
  end
end


local terra set_mappers()

end

--Single function for compiling the DSL
function compile_DSL( main_task, executable_name )

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
  local exe = os.getenv('OBJNAME') or executable_name
  if DSL_settings.mapper_path ~= nil then
    local cmapper, mapper_so = compile_mapper_run()
    print( "Mapper so file should be at: " .. mapper_so)
    regentlib.linklibrary(mapper_so)
    local link_flags = terralib.newlist({"-g", "-L" .. out_dir, "-lm", "-lhdf5", "-lmapper"})
    regentlib.saveobj(main_task, exe, "executable", cmapper.register_mappers, link_flags)
  else
   local link_flags = terralib.newlist({"-g", "-L" .. out_dir, "-lm", "-lhdf5"})
    regentlib.saveobj(main_task, exe, "executable", set_mappers, link_flags)
  end
end

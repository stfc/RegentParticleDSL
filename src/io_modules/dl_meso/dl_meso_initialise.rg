import "regent"

--require("src/io_modules/dl_meso/
require("src/RegentParticleDSL")

local zero_config = require("src/config/zero_config")
local zero_config_func = zero_config.generate_zero_config_quote(variables.config)
local zero_part_task = generate_zero_part_func()

local format = require("std/format")
local c = regentlib.c
local io_module = {}
function io_module.initialise(variables)


    local init_kernel = rquote
        var [variables.config] = region(ispace(int1d, 1), config_type)
        zero_config_func(variables.config);
        --End of initialisation (done by zero_config_func entirely now)
        dl_meso_read_mod.scan_control([variables.config])
        dl_meso_timing_mod.timchk([variables.config]) --TODO: NYI
        dl_meso_write_mod.write_output_header([variables.config])
        dl_meso_config_mod.sysdef(variables.config)
        var [variables.particle_array] = region(ispace(int1d, variables.config[0].nsyst), part)
        zero_part_task(variables.particle_array)
        for j in variables.particle_array.ispace do
            for i = 0,25 do
                [variables.particle_array][j].atmnam[i] = int8(0)
            end
        end
-- c.legion_runtime_issue_execution_fence(__runtime(), __context())
        dl_meso_config_mod.zero(variables.config, variables.particle_array)
-- c.legion_runtime_issue_execution_fence(__runtime(), __context())
        dl_meso_start_mod.start(variables.particle_array, variables.config)
-- c.legion_runtime_issue_execution_fence(__runtime(), __context())
        dl_meso_timing_mod.timchk([variables.config]) --TODO: NYI
-- c.legion_runtime_issue_execution_fence(__runtime(), __context())
        if [variables.config][0].kres == 0 and variables.config[0].straj == 0 and variables.config[0].ltraj then
            if variables.config[0].nseql > 0 then
                dl_meso_write_mod.gather_write_data(false, true, -variables.config[0].tstep * double(variables.config[0].nseql), variables.config, variables.particle_array)
            else
                dl_meso_write_mod.gather_write_data(false, true, 0.0, variables.config, variables.particle_array)
            end 
        end
    
        if variables.config[0].nstep == 0 then
            dl_meso_field_mod.plcfor_initial(variables.config, variables.particle_array)
            --IF (l_initial .AND. (.NOT. (l_config .AND. levcfg==2 .AND. nfold==1))) CALL write_config ('CFGINI', 2) --TODO NYI
-- c.legion_runtime_issue_execution_fence(__runtime(), __context())
            dl_meso_statistics_mod.statis(variables.config, variables.particle_array)
        end


        dl_meso_timing_mod.timchk([variables.config]) --TODO: NYI
        dl_meso_write_mod.write_output_summary(0.0, true, variables.config[0].l_scr, variables.config)
        
        --Set max cutoff
        for j in variables.particle_array.ispace do 
            variables.particle_array[j].core_part_space.cutoff = variables.config[0].rtcut
        end
    end
    return init_kernel
end


return io_module
--task main()
--    [initialise(variables)];
--end
--regentlib.start(main)

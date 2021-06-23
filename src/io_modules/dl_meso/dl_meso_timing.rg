import "regent"

local timing_mod = {}

--Timing functionality for dl_meso
--TODO: NYI
task timing_mod.timchk( config : region(ispace(int1d), config_type) ) where writes(config.tzero, config.time), reads(config.tzero, config.time)
do
    --TODO NYI
end

return timing_mod

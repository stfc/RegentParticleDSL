import "regent"

require("src/io_modules/dl_meso/dl_meso_config")
require("src/io_modules/dl_meso/dl_meso_particle")
require("src/io_modules/dl_meso/read_dl_meso")
require("src/io_modules/dl_meso/write_dl_meso")
local io_utils = require("src/io_modules/dl_meso/io_utils")
local format = require("std/format")
local c_unistd = terralib.includec("unistd.h")


task main()

--var cutoff : double = 0.0
--var srfzcut : double = 1.0
--var mxprm : int = 50
--var ltabpot : bool = true
--var lrcut : bool = true
--var nspe : int
--var namtmp : &&int8
--var masstmp : &double
--var chgetmp : &double
--var npot : int = 0
--var eunit : double = 1.0
--var gamma : &double = [&double](regentlib.c.malloc([terralib.sizeof(double)] * mxprm))

var i : int = 0
var j : int = 0

var config_region = region(ispace(int1d, 1), config_type)
config_region[0].cutoff = 1.0
config_region[0].mxprm = 0
config_region[0].ltabpot = false
config_region[0].lrcut = true
config_region[0].nspe = 0
config_region[0].npot = 0
for i = 0, CONST_max_species do
    for j = 0, 9 do
        config_region[0].namtmp[i][j] = int8(0)
    end
end
for i = 0, CONST_max_species do
    config_region[0].masstmp[i] = 0.0
    config_region[0].chgetmp[i] = 0.0
    config_region[0].ktype[i] = 0
end
config_region[0].eunit = 1.0
for i = 0, CONST_max_mxprm do
    config_region[0].gamma[i] = 0.0
    for j = 0, CONST_max_species * (CONST_max_species+1) / 2 do
        config_region[0].vvv[i][j] = 0.0
    end
end

--scan_field( &cutoff, srfzcut, &mxprm, ltabpot, lrcut, &nspe, &namtmp, &masstmp, &chgetmp )
scan_field(config_region)

regentlib.c.legion_runtime_issue_execution_fence(__runtime(), __context())
--c_unistd.sleep(5)
format.println("namtmp[0] is {}", [rawstring](config_region[0].namtmp[0]))
--
--read_field(nspe, masstmp, chgetmp, npot, mxprm, namtmp, eunit, gamma) 
read_field(config_region)

scan_control(config_region)
var a = 0
if config_region[0].l_exist then
    a = 1
end
format.println("l_exist is {}", a)

read_control(config_region, true, true, true)


--write_history_header(config_region)
end

regentlib.start(main)


import "regent"


local fkt = 2.0 / 3.0

local statistics_mod = {}

local sqrt = regentlib.sqrt(double)

local terra sclsum( sum_size : int, array : &double, stride : int) : double
    var k = 0
    var sum = 0.0
    for j=0, sum_size do
        sum = sum + array[k]
        k = k + stride
    end
    return sum
end

task statistics_mod.statis(config : region(ispace(int1d), config_type), parts : region(ispace(int1d), part)) where reads(config), writes(config), reads(parts), 
                                                                                        writes(parts.core_part_space)
do

--Unneeded: pass thermodynamic data to all nodes

    var rnf : double = 3.0 / (3.0 * config[0].nsyst) --TODO NYI: Constraints/frozen count 3.0_dp / (3.0_dp * REAL (nsyst-nfsyst, KIND=dp) - REAL(numcon, KIND=dp))
    config[0].pe = config[0].pe * rnf
    config[0].vir = config[0].vir * rnf
    config[0].ee = config[0].ee * rnf
    config[0].se = config[0].se * rnf
    --TODO NYI
    --If bonds or constraints
        --config[0].be = config[0].be * rnf
        --config[0].bdlng = config[0].bdlnd / double(config[0].numbond + config[0].numcon)
    --end
    --If angles
        --config[0].ae = config[0].ae * rnf
        --config[0].bdang = config[0].bdang / double(config[0].numang)
    --end
    --If dihedrals
        --config[0].de = config[0].de * rnf
        --config[0].bddhd = config[0].bddhd / double(config[0].numdhd)
    --end
    config[0].tke[0] = config[0].tke[0] * rnf
    config[0].tke[1] = config[0].tke[1] * rnf
    config[0].tke[2] = config[0].tke[2] * rnf
    for i = 0, 36 do
        config[0].stress[i] = config[0].stress[i] / config[0].volm
    end

    --final step values of all thermodynamic parameters
    config[0].stppe = config[0].pe --TODO NYI: Correction term + config[0].clr[0] / config[0].volm
    config[0].stpvir = config[0].vir --TODO NYI: Correction term + config[0].clr[1] / config[0].volm
    config[0].stptke = config[0].tke[0] + config[0].tke[1] + config[0].tke[2]
    config[0].stptkex = config[0].tke[0]
    config[0].stptkey = config[0].tke[1]
    config[0].stptkez = config[0].tke[2]
    config[0].stpte = config[0].stppe + config[0].stptke
    config[0].stpvlm = config[0].volm
    config[0].stpprs = double(config[0].nsyst) * (2.0 * config[0].stptke - config[0].stpvir) / (3.0 * config[0].stpvlm) --TODO NYI: nsyst - nfsyst
    config[0].stpttp = fkt * config[0].stptke
    config[0].stptpx = 2.0 * config[0].stptkex
    config[0].stptpy = 2.0 * config[0].stptkey
    config[0].stptpz = 2.0 * config[0].stptkez
    config[0].stpee = config[0].ee
    config[0].stpse = config[0].se
    config[0].stpbe = config[0].be
    config[0].stpbdl = config[0].bdlng
    config[0].stpbdmx = config[0].bdlmax
    config[0].stpbdmn = config[0].bdlmin
    config[0].stpae = config[0].ae
    config[0].stpang = config[0].bdang
    config[0].stpde = config[0].de
    config[0].stpdhd = config[0].bddhd
    var tempval = config[0].stress[32] + config[0].stress[33] + config[0].stress[34] + config[0].stress[35] 
                    - 0.5 * (config[0].stress[0] + config[0].stress[1] + config[0].stress[2] + config[0].stress[3]
                    + config[0].stress[16] + config[0].stress[17] + config[0].stress[18] + config[0].stress[19])
    config[0].stpzts = config[0].space.dim_z * tempval

    if config[0].nstep == 0 then
        --provide temporary values for rolling averages at timestep 0
        config[0].rav[0] = config[0].stpte
        config[0].rav[1] = config[0].stppe
        config[0].rav[2] = config[0].stpee
        config[0].rav[3] = config[0].stpse
        config[0].rav[4] = config[0].stpbe
        config[0].rav[5] = config[0].stpae
        config[0].rav[6] = config[0].stpde
        config[0].rav[7] = config[0].stpvir
        config[0].rav[8] = config[0].stptke
        config[0].rav[9] = config[0].stpprs
        config[0].rav[10] = config[0].stpvlm
        config[0].rav[11] = config[0].stpzts
        config[0].rav[12] = config[0].stpttp
        config[0].rav[13] = config[0].stptpx
        config[0].rav[14] = config[0].stptpy
        config[0].rav[15] = config[0].stptpz
    else
        --store quantities in stack
        var kstk = ((config[0].nstep-1) % config[0].nstk) + 1
        if config[0].nstep > config[0].nstk then
            if kstk == 1 then
                --Can't take a pointer to element of a region nor create variable sized array inputs so we cast as a workaround
                config[0].zum[0] = sclsum(config[0].nstk, [&double](config[0].stkpe), 1)
                config[0].zum[1] = sclsum(config[0].nstk, [&double](config[0].stkee), 1)
                config[0].zum[2] = sclsum(config[0].nstk, [&double](config[0].stkse), 1)
                config[0].zum[3] = sclsum(config[0].nstk, [&double](config[0].stkbe), 1)
                config[0].zum[4] = sclsum(config[0].nstk, [&double](config[0].stkae), 1)
                config[0].zum[5] = sclsum(config[0].nstk, [&double](config[0].stkde), 1)
                config[0].zum[6] = sclsum(config[0].nstk, [&double](config[0].stkvir), 1)
                config[0].zum[7] = sclsum(config[0].nstk, [&double](config[0].stkvlm), 1)
                config[0].zum[8] = sclsum(config[0].nstk, [&double](config[0].stkzts), 1)
                config[0].zum[9] = sclsum(config[0].nstk, [&double](config[0].stktkex), 1)
                config[0].zum[10] = sclsum(config[0].nstk, [&double](config[0].stktkey), 1)
                config[0].zum[11] = sclsum(config[0].nstk, [&double](config[0].stktkez), 1)
            end

            config[0].zum[0] = config[0].zum[0] - config[0].stkpe[kstk]
            config[0].zum[1] = config[0].zum[1] - config[0].stkee[kstk]
            config[0].zum[2] = config[0].zum[2] - config[0].stkse[kstk]
            config[0].zum[3] = config[0].zum[3] - config[0].stkbe[kstk]
            config[0].zum[4] = config[0].zum[4] - config[0].stkae[kstk]
            config[0].zum[5] = config[0].zum[5] - config[0].stkde[kstk]
            config[0].zum[6] = config[0].zum[6] - config[0].stkvir[kstk]
            config[0].zum[7] = config[0].zum[7] - config[0].stkvlm[kstk]
            config[0].zum[8] = config[0].zum[8] - config[0].stkzts[kstk]
            config[0].zum[9] = config[0].zum[9] - config[0].stktkex[kstk]
            config[0].zum[10] = config[0].zum[10] - config[0].stktkey[kstk]
            config[0].zum[11] = config[0].zum[11] - config[0].stktkez[kstk]

        end

        config[0].stkpe[kstk] = config[0].stppe
        config[0].stkvir[kstk] = config[0].stpvir
        config[0].stktkex[kstk] = config[0].stptkex
        config[0].stktkey[kstk] = config[0].stptkey
        config[0].stktkez[kstk] = config[0].stptkez
        config[0].stkee[kstk] = config[0].stpee
        config[0].stkse[kstk] = config[0].stpse
        config[0].stkbe[kstk] = config[0].stpbe
        config[0].stkae[kstk] = config[0].stpae
        config[0].stkde[kstk] = config[0].stpde
        config[0].stkvlm[kstk] = config[0].stpvlm
        config[0].stkzts[kstk] = config[0].stpzts
        config[0].zum[0] = config[0].zum[0] + config[0].stppe
        config[0].zum[1] = config[0].zum[1] + config[0].stpee
        config[0].zum[2] = config[0].zum[2] + config[0].stpse
        config[0].zum[3] = config[0].zum[3] + config[0].stpbe
        config[0].zum[4] = config[0].zum[4] + config[0].stpae
        config[0].zum[5] = config[0].zum[5] + config[0].stpde
        config[0].zum[6] = config[0].zum[6] + config[0].stpvir
        config[0].zum[7] = config[0].zum[7] + config[0].stpvlm
        config[0].zum[8] = config[0].zum[8] + config[0].stpzts
        config[0].zum[9] = config[0].zum[9] + config[0].stptkex
        config[0].zum[10] = config[0].zum[10] + config[0].stptkey
        config[0].zum[11] = config[0].zum[11] + config[0].stptkez

        --calculate rolling averages
        var zistk : double = double( regentlib.fmin(config[0].nstk, config[0].nstep))
        var rzistk : double = 1.0 / zistk
        var zumtke : double = config[0].zum[9] + config[0].zum[10] + config[0].zum[11]
        config[0].rav[0] = (config[0].zum[0] + zumtke) * rzistk
        config[0].rav[1] = config[0].zum[0] * rzistk
        config[0].rav[2] = config[0].zum[1] * rzistk
        config[0].rav[3] = config[0].zum[2] * rzistk
        config[0].rav[4] = config[0].zum[3] * rzistk
        config[0].rav[5] = config[0].zum[4] * rzistk
        config[0].rav[6] = config[0].zum[5] * rzistk
        config[0].rav[7] = config[0].zum[6] * rzistk
        config[0].rav[8] = zumtke * rzistk
        config[0].rav[9] = double(config[0].nsyst) * (2.0 * zumtke - config[0].zum[6]) / (3.0 * config[0].zum[7])
        config[0].rav[10] = config[0].zum[7] * rzistk
        config[0].rav[11] = config[0].zum[8] * rzistk
        config[0].rav[12] = fkt * zumtke * rzistk
        config[0].rav[13] = 2.0 * config[0].zum[9] * rzistk
        config[0].rav[14] = 2.0 * config[0].zum[10] * rzistk
        config[0].rav[15] = 2.0 * config[0].zum[11] * rzistk

        --Accumulate totals over steps
        config[0].nav = config[0].nav + 1
        var sclnv1 = double(config[0].nav - 1) / double(config[0].nav)
        var sclnv2 = 1.0 / double(config[0].nav)
        config[0].flc[0] = sclnv1 * (config[0].flc[0] + sclnv2 * (config[0].stpte - config[0].ave[0]) * (config[0].stpte - config[0].ave[0]))
        config[0].flc[1] = sclnv1 * (config[0].flc[1] + sclnv2 * (config[0].stppe - config[0].ave[1]) * (config[0].stppe - config[0].ave[1]))
        config[0].flc[2] = sclnv1 * (config[0].flc[2] + sclnv2 * (config[0].stpee - config[0].ave[2]) * (config[0].stpee - config[0].ave[2]))
        config[0].flc[3] = sclnv1 * (config[0].flc[3] + sclnv2 * (config[0].stpse - config[0].ave[3]) * (config[0].stpse - config[0].ave[3]))
        config[0].flc[4] = sclnv1 * (config[0].flc[4] + sclnv2 * (config[0].stpbe - config[0].ave[4]) * (config[0].stpbe - config[0].ave[4]))
        config[0].flc[5] = sclnv1 * (config[0].flc[5] + sclnv2 * (config[0].stpae - config[0].ave[5]) * (config[0].stpae - config[0].ave[5]))
        config[0].flc[6] = sclnv1 * (config[0].flc[6] + sclnv2 * (config[0].stpde - config[0].ave[6]) * (config[0].stpde - config[0].ave[6]))
        config[0].flc[7] = sclnv1 * (config[0].flc[7] + sclnv2 * (config[0].stpvir - config[0].ave[7]) * (config[0].stpvir - config[0].ave[7]))
        config[0].flc[8] = sclnv1 * (config[0].flc[8] + sclnv2 * (config[0].stptke - config[0].ave[8]) * (config[0].stptke - config[0].ave[8]))
        config[0].flc[9] = sclnv1 * (config[0].flc[9] + sclnv2 * (config[0].stpprs - config[0].ave[9]) * (config[0].stpprs - config[0].ave[9]))
        config[0].flc[10] = sclnv1 * (config[0].flc[10] + sclnv2 * (config[0].stpvlm - config[0].ave[10]) * (config[0].stpvlm - config[0].ave[10]))
        config[0].flc[11] = sclnv1 * (config[0].flc[11] + sclnv2 * (config[0].stpzts - config[0].ave[11]) * (config[0].stpzts - config[0].ave[11]))
        config[0].flc[12] = sclnv1 * (config[0].flc[12] + sclnv2 * (config[0].stpttp - config[0].ave[12]) * (config[0].stpttp - config[0].ave[12]))
        config[0].flc[13] = sclnv1 * (config[0].flc[13] + sclnv2 * (config[0].stptpx - config[0].ave[13]) * (config[0].stptpx - config[0].ave[13]))
        config[0].flc[14] = sclnv1 * (config[0].flc[14] + sclnv2 * (config[0].stptpy - config[0].ave[14]) * (config[0].stptpy - config[0].ave[14]))
        config[0].flc[15] = sclnv1 * (config[0].flc[15] + sclnv2 * (config[0].stptpz - config[0].ave[15]) * (config[0].stptpz - config[0].ave[15]))
        for i = 0, 36 do
            config[0].flc[16+i] = sclnv1 * (config[0].flc[16+i] + sclnv2 * (config[0].stress[i] - config[0].ave[16+i]) * 
                                    (config[0].stress[i] - config[0].ave[16+i]))
        end
        
        config[0].ave[0] = sclnv1 * config[0].ave[0] + sclnv2 * config[0].stpte
        config[0].ave[1] = sclnv1 * config[0].ave[1] + sclnv2 * config[0].stppe
        config[0].ave[2] = sclnv1 * config[0].ave[2] + sclnv2 * config[0].stpee
        config[0].ave[3] = sclnv1 * config[0].ave[3] + sclnv2 * config[0].stpse
        config[0].ave[4] = sclnv1 * config[0].ave[4] + sclnv2 * config[0].stpbe
        config[0].ave[5] = sclnv1 * config[0].ave[5] + sclnv2 * config[0].stpae
        config[0].ave[6] = sclnv1 * config[0].ave[6] + sclnv2 * config[0].stpde
        config[0].ave[7] = sclnv1 * config[0].ave[7] + sclnv2 * config[0].stpvir
        config[0].ave[8] = sclnv1 * config[0].ave[8] + sclnv2 * config[0].stptke
        config[0].ave[9] = sclnv1 * config[0].ave[9] + sclnv2 * config[0].stpprs
        config[0].ave[10] = sclnv1 * config[0].ave[10] + sclnv2 * config[0].stpvlm
        config[0].ave[11] = sclnv1 * config[0].ave[11] + sclnv2 * config[0].stpzts
        config[0].ave[12] = sclnv1 * config[0].ave[12] + sclnv2 * config[0].stpttp
        config[0].ave[13] = sclnv1 * config[0].ave[13] + sclnv2 * config[0].stptpx
        config[0].ave[14] = sclnv1 * config[0].ave[14] + sclnv2 * config[0].stptpy
        config[0].ave[15] = sclnv1 * config[0].ave[15] + sclnv2 * config[0].stptpz
        for i = 0, 36 do
            config[0].ave[16+i] = sclnv1* config[0].ave[16+i] + sclnv2 * config[0].stress[i]
        end

    end

    --temperature scaling
    if config[0].nstep > config[0].nseql then
        config[0].ltemp = false
    end

    if config[0].nstep <= config[0].nseql then
        if config[0].ltemp then
            if config[0].nstep % config[0].nsbts == 0 then
                var tscal = sqrt(config[0].temp / config[0].stpttp)
                for i in parts.ispace do
                    parts[i].core_part_space.vel_x = tscal * parts[i].core_part_space.vel_x
                    parts[i].core_part_space.vel_y = tscal * parts[i].core_part_space.vel_y
                    parts[i].core_part_space.vel_z = tscal * parts[i].core_part_space.vel_z
                end
                --TODO NYI: Constraints and quenching
            end
        end
        config[0].nav = 0
        for i = 0, CONST_statsize do
            config[0].ave[i] = 0.0
            config[0].flc[i] = 0.0
        end
    end

end

return statistics_mod

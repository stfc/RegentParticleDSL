import "regent"

local rt12=3.464101615377546 --!< Square root of 12 used for approximation of Gaussian random

local c = regentlib.c

local field_mod = {}

local function force_mdvv(part1, part2, r2, config)
    local sqrt = regentlib.sqrt(double)
    local kernel = rquote
        var xdif = part1.core_part_space.pos_x - part2.core_part_space.pos_x
        var ydif = part1.core_part_space.pos_y - part2.core_part_space.pos_y
        var zdif = part1.core_part_space.pos_z - part2.core_part_space.pos_z
        var rrr = sqrt(r2)
        var ai = double(part1.ltp)
        var aj = double(part2.ltp)
        var ab : double
        if ai > aj then
            ab = ai * (ai -1.0) * 0.5 + aj + 0.5
        else
            ab = aj * (aj - 1.0) * 0.5 + ai + 0.5
        end
        var k : int = int(ab) --TODO: This may be wrong - needs checking since we use 0-index instead of 1 index for both ltp and k

        --calculate conservative interaction force and potential energy if TODO: NYI not excluding interactions in bond,constraint or angle
        var gforce = 0.0
        var pot = 0.0
        --If not excluded then
            if config.ktype[k] == 2 then
                --espanol and warren potential
                if r2 < config.vvv[2][k] then
                    var vk0 : double = config.vvv[0][k]
                    var vk1 : double = config.vvv[1][k]
                    var vk2 : double = 1.0 / vk1
                    var scrn : double = 1.0 / rrr - vk2
                    pot = 0.5 * vk0 * (vk1 - rrr) * (vk1 - rrr) * vk2
                    gforce = vk0 * scrn
                end
            else
                --Other ktypes not yet supported -- remove this assert for performance runs
                regentlib.assert(false, "Not supporting other ktypes yet")
            end
            --Accumulate potential energy and surface energy
            config.pe += pot
            --TODO NYI: config[0].se += config[0].lsurf[k] * pot            
        --end
        --calculate random and drag forces (dpd thermostat)
        var rforce : double = 0.0
        var dforce : double = 0.0
        if r2 < config.rtct2 then
            var vxdif : double = part1.core_part_space.vel_x - part2.core_part_space.vel_x
            var vydif : double = part1.core_part_space.vel_y - part2.core_part_space.vel_y
            var vzdif : double = part1.core_part_space.vel_z - part2.core_part_space.vel_z
            var ran : double = rt12 * (c.drand48() - 0.5)  --FIXME: drand48 inside parallel code is not ideal
            var scrn : double = 1.0 / rrr - config.rrtcut
            rforce = config.sigma[k] * scrn * ran
            dforce = -config.gamma[k] * scrn * scrn * (xdif * vxdif + ydif * vydif + zdif * vzdif)
        end
        --sum of forces
        var sumforce : double = gforce + rforce + dforce
        -- assign virial term
        config.vir -= sumforce * r2

        --assign stress terms
        config.mdvv_type.strscxx += gforce * xdif * xdif
        config.mdvv_type.strscxy += gforce * xdif * ydif
        config.mdvv_type.strscxz += gforce * xdif * zdif
        config.mdvv_type.strscyy += gforce * ydif * ydif
        config.mdvv_type.strscyz += gforce * ydif * zdif
        config.mdvv_type.strsczz += gforce * zdif * zdif

        config.mdvv_type.strsdxx += dforce * xdif * xdif
        config.mdvv_type.strsdxy += dforce * xdif * ydif
        config.mdvv_type.strsdxz += dforce * xdif * zdif
        config.mdvv_type.strsdyy += dforce * ydif * ydif
        config.mdvv_type.strsdyz += dforce * ydif * zdif
        config.mdvv_type.strsdzz += dforce * zdif * zdif

        config.mdvv_type.strsrxx += rforce * xdif * xdif
        config.mdvv_type.strsrxy += rforce * xdif * ydif
        config.mdvv_type.strsrxz += rforce * xdif * zdif
        config.mdvv_type.strsryy += rforce * ydif * ydif
        config.mdvv_type.strsryz += rforce * ydif * zdif
        config.mdvv_type.strsrzz += rforce * zdif * zdif

        part1.fxx -= sumforce * xdif
        part1.fyy -= sumforce * ydif
        part1.fzz -= sumforce * zdif
        part2.fxx += sumforce * xdif
        part2.fyy += sumforce * ydif
        part2.fzz += sumforce * zdif

    end

    return kernel

end

task field_mod.forces_mdvv(config : region(ispace(int1d), config_type), parts : region(ispace(int1d), part)) where writes(config), reads(config), writes(parts), reads(parts)
do

    var ab : double = 0.0
    var ai : double = 0.0
    var aj : double = 0.0
    var gforce : double = 0.0
    var pot : double = 0.0
    var rrr : double = 0.0
    var rsq : double = 0.0
    var ran : double = 0.0
    var rvolm : double = 0.0
    var xdif : double = 0.0
    var ydif : double = 0.0
    var zdif : double = 0.0
    var vxdif : double = 0.0
    var vydif : double = 0.0
    var vzdif : double = 0.0
    var rforce : double = 0.0
    var dforce : double = 0.0
    var scrn : double = 0.0
    var sumforce : double = 0.0
    config[0].mdvv_type.strscxx = 0.0
    config[0].mdvv_type.strscxy = 0.0
    config[0].mdvv_type.strscxz = 0.0
    config[0].mdvv_type.strscyy = 0.0
    config[0].mdvv_type.strscyz = 0.0
    config[0].mdvv_type.strsczz = 0.0
    config[0].mdvv_type.strsdxx = 0.0
    config[0].mdvv_type.strsdxy = 0.0
    config[0].mdvv_type.strsdxz = 0.0
    config[0].mdvv_type.strsdyy = 0.0
    config[0].mdvv_type.strsdyz = 0.0
    config[0].mdvv_type.strsdzz = 0.0
    config[0].mdvv_type.strsrxx = 0.0
    config[0].mdvv_type.strsrxy = 0.0
    config[0].mdvv_type.strsrxz = 0.0
    config[0].mdvv_type.strsryy = 0.0
    config[0].mdvv_type.strsryz = 0.0
    config[0].mdvv_type.strsrzz = 0.0

    for i in parts.ispace do
        parts[i].fxx = 0.0
        parts[i].fyy = 0.0
        parts[i].fzz = 0.0
    end

    --[invoke(rexpr config[0] end, {force_mdvv, SYMMETRIC_PAIRWISE}, BARRIER)];
    --TODO: Would be much better to do this in parallel but cannot use invoke here at this time.
    for i in parts.ispace do
        for j in parts.ispace do
            if i < j then
                var dx = parts[i].core_part_space.pos_x - parts[j].core_part_space.pos_x
                var dy = parts[i].core_part_space.pos_y - parts[j].core_part_space.pos_y
                var dz = parts[i].core_part_space.pos_z - parts[j].core_part_space.pos_z
                var r2 = dx*dx + dy*dy + dz*dz
                if r2 <= (config[0].cutoff*config[0].cutoff) then
                    [force_mdvv(rexpr parts[i] end, rexpr parts[j] end, rexpr r2 end, rexpr config[0] end)]
                end
            end
        end
    end
    
    rvolm = 1.0 / (3.0 * config[0].volm) --Multinode (3.0_dp * volm * REAL(nodes, KIND=dp))
    
    config[0].ivrl[0] = config[0].ivrl[0] - config[0].mdvv_type.strscxx --TODO NYI: LJ correction + clr(2) * rvolm
    config[0].ivrl[1] = config[0].ivrl[1] - config[0].mdvv_type.strscyy --TODO NYI: LJ correction + clr(2) * rvolm
    config[0].ivrl[2] = config[0].ivrl[2] - config[0].mdvv_type.strsczz --TODO NYI: LJ correction + clr(2) * rvolm

    rvolm = 1.0 / 3.0 --Multinode (3.0_dp * REAL(nodes, KIND=dp)
    config[0].stress[0] = config[0].stress[0] + config[0].mdvv_type.strscxx --TODO NYI: LJ correction - clr(2) * rvolm
    config[0].stress[1] = config[0].stress[1] + config[0].mdvv_type.strsdxx
    config[0].stress[2] = config[0].stress[2] + config[0].mdvv_type.strsrxx

    config[0].stress[4] = config[0].stress[4] + config[0].mdvv_type.strscxy
    config[0].stress[5] = config[0].stress[5] + config[0].mdvv_type.strsdxy
    config[0].stress[6] = config[0].stress[6] + config[0].mdvv_type.strsrxy
    
    config[0].stress[8] = config[0].stress[8] + config[0].mdvv_type.strscxz
    config[0].stress[9] = config[0].stress[9] + config[0].mdvv_type.strsdxz
    config[0].stress[10] = config[0].stress[10] + config[0].mdvv_type.strsrxz

    config[0].stress[12] = config[0].stress[12] + config[0].mdvv_type.strscxy
    config[0].stress[13] = config[0].stress[13] + config[0].mdvv_type.strsdxy
    config[0].stress[14] = config[0].stress[14] + config[0].mdvv_type.strsrxy

    config[0].stress[16] = config[0].stress[16] + config[0].mdvv_type.strscyy --TODO NYI: LJ correction - clr(2) * rvolm
    config[0].stress[17] = config[0].stress[17] + config[0].mdvv_type.strsdyy
    config[0].stress[18] = config[0].stress[18] + config[0].mdvv_type.strsryy

    config[0].stress[20] = config[0].stress[20] + config[0].mdvv_type.strscyz
    config[0].stress[21] = config[0].stress[21] + config[0].mdvv_type.strsdyz
    config[0].stress[22] = config[0].stress[22] + config[0].mdvv_type.strsryz

    config[0].stress[24] = config[0].stress[24] + config[0].mdvv_type.strscxz
    config[0].stress[25] = config[0].stress[25] + config[0].mdvv_type.strsdxz
    config[0].stress[26] = config[0].stress[26] + config[0].mdvv_type.strsrxz

    config[0].stress[28] = config[0].stress[28] + config[0].mdvv_type.strscyz
    config[0].stress[29] = config[0].stress[29] + config[0].mdvv_type.strsdyz
    config[0].stress[30] = config[0].stress[30] + config[0].mdvv_type.strsryz

    config[0].stress[32] = config[0].stress[32] + config[0].mdvv_type.strsczz --TODO NYI: LJ correction - clr(2) * rvolm
    config[0].stress[33] = config[0].stress[33] + config[0].mdvv_type.strsdzz 
    config[0].stress[34] = config[0].stress[34] + config[0].mdvv_type.strsrzz 
end



task field_mod.plcfor_initial(config : region(ispace(int1d), config_type), parts : region(ispace(int1d), part)) where writes(config), reads(config), writes(parts), reads(parts)
do

    config[0].pe = 0.0
    config[0].ee = 0.0
    config[0].se = 0.0
    config[0].be = 0.0
    config[0].ae = 0.0
    config[0].de = 0.0
    config[0].vir = 0.0
    config[0].ivrl[0] = 0.0
    config[0].ivrl[1] = 0.0
    config[0].ivrl[2] = 0.0
    for i = 0, 36 do
        config[0].stress[i] = 0.0
    end
    config[0].bdlng = 0.0
    config[0].bdlmin = config[0].volm
    config[0].bdlmax = 0.0
    config[0].bdang = 0.0
    config[0].bddhd = 0.0
    var strsxx : double = 0.0
    var strsxy : double = 0.0
    var strsxz : double = 0.0
    var strsyy : double = 0.0
    var strsyz : double = 0.0
    var strszz : double = 0.0

    --No need for loc/lmp/lblclst/nlist

    --TODO NYI: calculate localized densities for many-body dpd
    --TODO NYI : export particle data - likely unneeded
    --Unneeded: sort global/local bead list
    --Unneeded: construct neighbour list
    --determine if forces need to be calculate
    var l_force : bool = true --TODO: NYI ((.NOT. l_config) .OR. (l_config .AND. levcfg<2))
    if l_force then
        regentlib.assert(config[0].itype == 0, "itypes other than 0 are not yet supported")
        field_mod.forces_mdvv(config, parts)
    else
        regentlib.assert(false, "l_force is required")
    end

    --TODO NYI: determine potential energies for many-body dpd (line 3473)

    --TODO NYI: determine electrostatic forces lines 3477-3503

    --TODO NYI: determine all bonded forces (lines 3507-3521)

    --TODO NYI: determine surface/wall forces, potential and virial

    --Don't need to sum forces from halo regions.

    --TODO NYI: quench forces and velocities for frozen particles (lines 3538-3540)

    --Calculate kinetic energy (divided into components) and stress tensor
    strsxx = 0.0 
    strsxy = 0.0 
    strsxz = 0.0 
    strsyy = 0.0 
    strsyz = 0.0 
    strszz = 0.0 
    for i in parts.ispace do
        strsxx = strsxx + parts[i].core_part_space.mass * parts[i].core_part_space.vel_x * parts[i].core_part_space.vel_x
        strsxy = strsxy + parts[i].core_part_space.mass * parts[i].core_part_space.vel_x * parts[i].core_part_space.vel_y
        strsxz = strsxz + parts[i].core_part_space.mass * parts[i].core_part_space.vel_x * parts[i].core_part_space.vel_z
        strsyy = strsyy + parts[i].core_part_space.mass * parts[i].core_part_space.vel_y * parts[i].core_part_space.vel_y
        strsyz = strsyz + parts[i].core_part_space.mass * parts[i].core_part_space.vel_y * parts[i].core_part_space.vel_z
        strszz = strszz + parts[i].core_part_space.mass * parts[i].core_part_space.vel_z * parts[i].core_part_space.vel_z
    end
    config[0].tke[0] = strsxx * 0.5
    config[0].tke[1] = strsyy * 0.5
    config[0].tke[2] = strszz * 0.5

    config[0].stress[3] = strsxx
    config[0].stress[7] = strsxy
    config[0].stress[11] = strsxz
    config[0].stress[15] = strsxy
    config[0].stress[19] = strsyy
    config[0].stress[23] = strsyz
    config[0].stress[27] = strsxz
    config[0].stress[31] = strsyz
    config[0].stress[35] = strszz

end


return field_mod

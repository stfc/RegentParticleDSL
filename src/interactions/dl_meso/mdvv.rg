import "regent"

require("src/RegentParticleDSL")
import_dl_meso()

function reassign_particle_properties(part1, config)
    local kernel = rquote
        part1.core_part_space.mass = config.masstmp[part1.ltp]
--        for x = 0, 9 do
--            part1.atmnam[x] = config.namspe[part1.ltp][x]
--        end
    end
    return kernel
end

function deportdata(part1, config)
    local floor = regentlib.floor(double)
    local kernel = rquote
        --TODO NYI: apply lees-edwards shearing periodic boundaries
        part1.core_part_space.pos_x = part1.core_part_space.pos_x - config.space.dim_x *  floor(part1.core_part_space.pos_x / config.space.dim_x)
        part1.core_part_space.pos_y = part1.core_part_space.pos_y - config.space.dim_y *  floor(part1.core_part_space.pos_y / config.space.dim_y)
        part1.core_part_space.pos_z = part1.core_part_space.pos_z - config.space.dim_z *  floor(part1.core_part_space.pos_z / config.space.dim_z)
    end
    return kernel
end

function velocity_verlet_stage1(part1, config)
    local floor = regentlib.floor(double)
    local kernel = rquote
        --TODO NYI: Only non-frozen beads
        --var qi =  config.chgetmp[part1.ltp] TODO NYI: Electrostatic field
        var rmass = 1.0 / part1.core_part_space.mass
        part1.core_part_space.vel_x = part1.core_part_space.vel_x + 0.5 * config.tstep * ((part1.fxx) * rmass) --TODO NYI: constant force bdfrcx
        part1.core_part_space.vel_y = part1.core_part_space.vel_y + 0.5 * config.tstep * ((part1.fyy) * rmass) --TODO NYI: constant force bdfrcx
        part1.core_part_space.vel_z = part1.core_part_space.vel_z + 0.5 * config.tstep * ((part1.fzz) * rmass) --TODO NYI: constant force bdfrcx
    
        part1.core_part_space.pos_x = part1.core_part_space.pos_x + config.tstep * part1.core_part_space.vel_x
        part1.core_part_space.pos_y = part1.core_part_space.pos_y + config.tstep * part1.core_part_space.vel_y
        part1.core_part_space.pos_z = part1.core_part_space.pos_z + config.tstep * part1.core_part_space.vel_z;

        --For now we put deportdata in here too
        [deportdata(part1, config)];

    end
    return kernel
end

function velocity_verlet_stage2(part1, config)
    local kernel = rquote
        --TODO NYI: Only non-frozen beads
        --var qi = config.chgetmp[part1.ltp] TODO NYI: Electrostatic field
        var rmass = 1.0 / part1.core_part_space.mass
        part1.core_part_space.vel_x = part1.core_part_space.vel_x + 0.5 * config.tstep * ((part1.fxx) * rmass) --TODO NYI: constant force bdfrcx
        part1.core_part_space.vel_y = part1.core_part_space.vel_y + 0.5 * config.tstep * ((part1.fyy) * rmass) --TODO NYI: constant force bdfrcx
        part1.core_part_space.vel_z = part1.core_part_space.vel_z + 0.5 * config.tstep * ((part1.fzz) * rmass) --TODO NYI: constant force bdfrcx
    end
    return kernel
end

function mdvv_kinetic_energy_compute(part1, config)
    local kernel = rquote
        config.mdvv_type.strscxx += part1.core_part_space.mass * part1.core_part_space.vel_x * part1.core_part_space.vel_x
        config.mdvv_type.strscxy += part1.core_part_space.mass * part1.core_part_space.vel_x * part1.core_part_space.vel_y
        config.mdvv_type.strscxz += part1.core_part_space.mass * part1.core_part_space.vel_x * part1.core_part_space.vel_z
        config.mdvv_type.strscyy += part1.core_part_space.mass * part1.core_part_space.vel_y * part1.core_part_space.vel_y
        config.mdvv_type.strscyz += part1.core_part_space.mass * part1.core_part_space.vel_y * part1.core_part_space.vel_z
        config.mdvv_type.strsczz += part1.core_part_space.mass * part1.core_part_space.vel_z * part1.core_part_space.vel_z
    end
    return kernel
end


function mdvv_nvt(stage)

    local mdvv_nvt_kern = rquote
        --TODO NYI: if using constraints, work out lists of beads to work on

        var strsxx : double = 0.0
        var strsxy : double = 0.0
        var strsxz : double = 0.0
        var strsyy : double = 0.0
        var strsyz : double = 0.0
        var strszz : double = 0.0

        if stage == 1 then
            --Move particles by velocity verlet algorithm
            --At the moment deportdata is sync'd in here
            [invoke(variables.config, {velocity_verlet_stage1, PER_PART}, NO_BARRIER)];
       
            --TODO NYI: Apply boundary conditions covered by srftype

            --TODO NYI: If using constraints, apply stage 1 of RATTLE 

            --TODO NYI: If using shear, apply sliding distance for shearing
            

        elseif stage == 2 then
            --Move particles by velocity verlet algorithm (this time velocity updates only)
            [invoke(variables.config, {velocity_verlet_stage2, PER_PART}, BARRIER)];

            --TODO NYI: If using constraints, apply stage 2 of RATTLE
    
            --calculate kinetic energy and kinetic part of stress tensor
            variables.config[0].mdvv_type.strscxx = 0.0
            variables.config[0].mdvv_type.strscxy = 0.0
            variables.config[0].mdvv_type.strscxz = 0.0
            variables.config[0].mdvv_type.strscyy = 0.0
            variables.config[0].mdvv_type.strscyz = 0.0
            variables.config[0].mdvv_type.strsczz = 0.0
            [invoke(variables.config, {mdvv_kinetic_energy_compute, PER_PART}, BARRIER)];
            variables.config[0].tke[0] = 0.5 * variables.config[0].mdvv_type.strscxx
            variables.config[0].tke[1] = 0.5 * variables.config[0].mdvv_type.strscyy
            variables.config[0].tke[2] = 0.5 * variables.config[0].mdvv_type.strsczz

            --assign stress components to full tensor
            variables.config[0].stress[3] = variables.config[0].mdvv_type.strscxx
            variables.config[0].stress[7] = variables.config[0].mdvv_type.strscxy
            variables.config[0].stress[11] = variables.config[0].mdvv_type.strscxz
            variables.config[0].stress[15] = variables.config[0].mdvv_type.strscxy
            variables.config[0].stress[19] = variables.config[0].mdvv_type.strscyy
            variables.config[0].stress[23] = variables.config[0].mdvv_type.strscyz
            variables.config[0].stress[27] = variables.config[0].mdvv_type.strscxz
            variables.config[0].stress[31] = variables.config[0].mdvv_type.strscyz
            variables.config[0].stress[35] = variables.config[0].mdvv_type.strsczz

        end        
    end

    return mdvv_nvt_kern
end


local sqrt = regentlib.sqrt(double)

function zero_forces(part1, config)
    local zero_force_kernel = rquote
        part1.fxx = 0.0
        part1.fyy = 0.0
        part1.fzz = 0.0
    end
    return zero_force_kernel
end

function forces_mdvv(part1, part2, r2, config)
    local rt12=3.464101615377546 --!< Square root of 12 used for approximation of Gaussian random
    local c = regentlib.c
    local floor = regentlib.floor(double)
    local forces_kernel = rquote
        var xdif : double = part2.core_part_space.pos_x - part1.core_part_space.pos_x
        var ydif : double = part2.core_part_space.pos_y - part1.core_part_space.pos_y
        var zdif : double = part2.core_part_space.pos_z - part1.core_part_space.pos_z
        if xdif > 0.5 * config.space.dim_x then
            xdif = xdif - config.space.dim_x
        end
        if xdif < -0.5 * config.space.dim_x then
            xdif = xdif + config.space.dim_x
        end
        if ydif > 0.5 * config.space.dim_y then
            ydif = ydif - config.space.dim_y
        end
        if ydif < -0.5 * config.space.dim_y then
            ydif = ydif + config.space.dim_y
        end
        if zdif > 0.5 * config.space.dim_z then
            zdif = zdif - config.space.dim_z
        end
        if zdif < -0.5 * config.space.dim_z then
            zdif = zdif + config.space.dim_z
        end
        var rrr : double = sqrt(r2)
        var ai : double = double(part1.ltp)
        var aj : double = double(part2.ltp)
        var ab : double
        if ai > aj then
            ab = ai * (ai -1.0) * 0.5 + aj + 0.5
        else
            ab = aj * (aj - 1.0) * 0.5 + ai + 0.5
        end
        var k : int = int(ab) --TODO: This may be wrong - needs checking since we use 0-index instead of 1 index for both ltp and k

        --Conservative force
        var gforce : double = 0.0
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
            end
            config.pe += pot
            --TODO NYI: config[0].se += config[0].lsurf[k] * pot
        --end
   
        --calculate random and drag forces (dpd thermostat)
        var rforce : double = 0.0
        var dforce : double = 0.0
        var scrn : double = 0.0
        if r2 < config.rtct2 then
            var vxdif : double = part2.core_part_space.vel_x - part1.core_part_space.vel_x
            var vydif : double = part2.core_part_space.vel_y - part1.core_part_space.vel_y
            var vzdif : double = part2.core_part_space.vel_z - part1.core_part_space.vel_z
            var ran : double = rt12 * (c.drand48() - 0.5)  --FIXME: drand48 inside parallel code is not ideal
            scrn = 1.0 / rrr - config.rrtcut
            rforce = config.sigma[k] * scrn * ran
            dforce = -config.gamma[k] * scrn * scrn * (xdif * vxdif + ydif * vydif + zdif * vzdif)
        end
        --sum of forces
        var sumforce : double = gforce + rforce + dforce 

        config.vir += -sumforce * r2
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

        -- assign force terms
--        if part1.lab == 1 then
--            format.println("{} {} fxx = {} {}",part1.lab, part2.lab, -sumforce * xdif, xdif)
--        end
        part1.fxx += -sumforce * xdif
        part1.fyy += -sumforce * ydif
        part1.fzz += -sumforce * zdif
        part2.fxx += sumforce * xdif
        part2.fyy += sumforce * ydif
        part2.fzz += sumforce * zdif        
    end
    return forces_kernel
end

function plcfor_mdvv()
    local plcfor_mdvv_kern = rquote
        --Reset accumulators
        variables.config[0].mdvv_type.strscxx = 0.0
        variables.config[0].mdvv_type.strscxy = 0.0
        variables.config[0].mdvv_type.strscxz = 0.0
        variables.config[0].mdvv_type.strscyy = 0.0
        variables.config[0].mdvv_type.strscyz = 0.0
        variables.config[0].mdvv_type.strsczz = 0.0
        variables.config[0].mdvv_type.strsdxx = 0.0
        variables.config[0].mdvv_type.strsdxy = 0.0
        variables.config[0].mdvv_type.strsdxz = 0.0
        variables.config[0].mdvv_type.strsdyy = 0.0
        variables.config[0].mdvv_type.strsdyz = 0.0
        variables.config[0].mdvv_type.strsdzz = 0.0
        variables.config[0].mdvv_type.strsrxx = 0.0
        variables.config[0].mdvv_type.strsrxy = 0.0
        variables.config[0].mdvv_type.strsrxz = 0.0
        variables.config[0].mdvv_type.strsryy = 0.0
        variables.config[0].mdvv_type.strsryz = 0.0
        variables.config[0].mdvv_type.strsrzz = 0.0;
        [invoke(variables.config, {zero_forces, PER_PART}, NO_BARRIER)];
        --Zero force values
        --TODO: NYI calculate localized densities for many-body dpd

        --Calculate pair forces
        --forces_mdvv
        [invoke(variables.config, {forces_mdvv, SYMMETRIC_PAIRWISE}, BARRIER)];
--        var sumx = 0.0
--        var sumy = 0.0
--        var sumz = 0.0
--        var velx = 0.0
--        var vely = 0.0
--        var velz = 0.0
--        for i in neighbour_init.padded_particle_array.ispace do
--            if neighbour_init.padded_particle_array[i].neighbour_part_space._valid then
--                sumx = sumx + neighbour_init.padded_particle_array[i].fxx
--                sumy = sumy + neighbour_init.padded_particle_array[i].fyy
--                sumz = sumz + neighbour_init.padded_particle_array[i].fzz
--                velx = neighbour_init.padded_particle_array[i].core_part_space.vel_x
--                vely = neighbour_init.padded_particle_array[i].core_part_space.vel_y
--                velz = neighbour_init.padded_particle_array[i].core_part_space.vel_z
--            end
--        end
--        format.println("Force sums are {} {} {}", sumx, sumy, sumz)
--        format.println("Velocity sums are {} {} {}", velx, vely, velz)
    	var rvolm = 1.0 / (3.0 * variables.config[0].volm) --Multinode (3.0_dp * volm * REAL(nodes, KIND=dp))

    	variables.config[0].ivrl[0] = variables.config[0].ivrl[0] - variables.config[0].mdvv_type.strscxx --TODO NYI: LJ correction + clr(2) * rvolm
    	variables.config[0].ivrl[1] = variables.config[0].ivrl[1] - variables.config[0].mdvv_type.strscyy --TODO NYI: LJ correction + clr(2) * rvolm
    	variables.config[0].ivrl[2] = variables.config[0].ivrl[2] - variables.config[0].mdvv_type.strsczz --TODO NYI: LJ correction + clr(2) * rvolm

    	rvolm = 1.0 / 3.0 --Multinode (3.0_dp * REAL(nodes, KIND=dp)
    	variables.config[0].stress[0] = variables.config[0].stress[0] + variables.config[0].mdvv_type.strscxx --TODO NYI: LJ correction - clr(2) * rvolm
    	variables.config[0].stress[1] = variables.config[0].stress[1] + variables.config[0].mdvv_type.strsdxx
    	variables.config[0].stress[2] = variables.config[0].stress[2] + variables.config[0].mdvv_type.strsrxx

    	variables.config[0].stress[4] = variables.config[0].stress[4] + variables.config[0].mdvv_type.strscxy
    	variables.config[0].stress[5] = variables.config[0].stress[5] + variables.config[0].mdvv_type.strsdxy
    	variables.config[0].stress[6] = variables.config[0].stress[6] + variables.config[0].mdvv_type.strsrxy

    	variables.config[0].stress[8] = variables.config[0].stress[8] + variables.config[0].mdvv_type.strscxz
    	variables.config[0].stress[9] = variables.config[0].stress[9] + variables.config[0].mdvv_type.strsdxz
    	variables.config[0].stress[10] = variables.config[0].stress[10] + variables.config[0].mdvv_type.strsrxz

    	variables.config[0].stress[12] = variables.config[0].stress[12] + variables.config[0].mdvv_type.strscxy
    	variables.config[0].stress[13] = variables.config[0].stress[13] + variables.config[0].mdvv_type.strsdxy
    	variables.config[0].stress[14] = variables.config[0].stress[14] + variables.config[0].mdvv_type.strsrxy

    	variables.config[0].stress[16] = variables.config[0].stress[16] + variables.config[0].mdvv_type.strscyy --TODO NYI: LJ correction - clr(2) * rvolm
    	variables.config[0].stress[17] = variables.config[0].stress[17] + variables.config[0].mdvv_type.strsdyy
    	variables.config[0].stress[18] = variables.config[0].stress[18] + variables.config[0].mdvv_type.strsryy

    	variables.config[0].stress[20] = variables.config[0].stress[20] + variables.config[0].mdvv_type.strscyz
    	variables.config[0].stress[21] = variables.config[0].stress[21] + variables.config[0].mdvv_type.strsdyz
    	variables.config[0].stress[22] = variables.config[0].stress[22] + variables.config[0].mdvv_type.strsryz

    	variables.config[0].stress[24] = variables.config[0].stress[24] + variables.config[0].mdvv_type.strscxz
    	variables.config[0].stress[25] = variables.config[0].stress[25] + variables.config[0].mdvv_type.strsdxz
    	variables.config[0].stress[26] = variables.config[0].stress[26] + variables.config[0].mdvv_type.strsrxz

    	variables.config[0].stress[28] = variables.config[0].stress[28] + variables.config[0].mdvv_type.strscyz
    	variables.config[0].stress[29] = variables.config[0].stress[29] + variables.config[0].mdvv_type.strsdyz
    	variables.config[0].stress[30] = variables.config[0].stress[30] + variables.config[0].mdvv_type.strsryz

    	variables.config[0].stress[32] = variables.config[0].stress[32] + variables.config[0].mdvv_type.strsczz --TODO NYI: LJ correction - clr(2) * rvolm
    	variables.config[0].stress[33] = variables.config[0].stress[33] + variables.config[0].mdvv_type.strsdzz
    	variables.config[0].stress[34] = variables.config[0].stress[34] + variables.config[0].mdvv_type.strsrzz

        --TODO: NYI determine potential energies for many-body dpd

        --TODO: NYI determine electrostatic forces

        --TODO: NYI determine all bonded forces

        --TODO: NYI determine surface/wall forces, potential and virial        

        --TODO: NYI sum forces and virials from halo region -- check if this has stuff we need, don't think so

        --TODO: NYI quench forces and velocities for frozen particles

    end
    return plcfor_mdvv_kern
end

local c_stdio = terralib.includec("stdio.h")

task mdvv()


[dl_meso_init.initialise(variables)];

[neighbour_init.initialise(variables)];
[neighbour_init.update_cells(variables)];

var finish : bool = false
var klock = 0
var l_hist : bool
var l_exp : bool


regentlib.assert(variables.config[0].btype == 0, "Only btype 0 supported")

format.println("Starting DL_MESO for {} steps", variables.config[0].nrun);
while not finish do

    variables.config[0].nstep = variables.config[0].nstep + 1
    --dl_meso_timing_mod.timchk(config) --TODO Include this.

    --Reset accumulators
    variables.config[0].pe = 0.0
    variables.config[0].ee = 0.0
    variables.config[0].se = 0.0
    variables.config[0].be = 0.0
    variables.config[0].ae = 0.0
    variables.config[0].de = 0.0
    variables.config[0].vir = 0.0
    variables.config[0].ivrl[0] = 0.0
    variables.config[0].ivrl[1] = 0.0
    variables.config[0].ivrl[2] = 0.0
    for i=0, 36 do
        variables.config[0].stress[i] = 0.0
    end
    variables.config[0].bdlng = 0.0
    variables.config[0].bdlmin = variables.config[0].volm
    variables.config[0].bdlmax = 0.0
    variables.config[0].bdang = 0.0
    variables.config[0].bddhd = 0.0

    if variables.config[0].btype == 0 then
        [mdvv_nvt(1)];
    end

    [invoke(variables.config, {reassign_particle_properties, PER_PART}, NO_BARRIER)];
    --dl_meso_timing_mod.timchk(config) --TODO Include this.
--    format.println("Step {}", klock);
    [plcfor_mdvv()]; 
--    for i in neighbour_init.padded_particle_array.ispace do
--        if neighbour_init.padded_particle_array[i].lab == 1 and neighbour_init.padded_particle_array[i].neighbour_part_space._valid then
--            format.println("part forces {} {} {} {}", neighbour_init.padded_particle_array[i].lab, 
--                                               neighbour_init.padded_particle_array[i].fxx,
--                                               neighbour_init.padded_particle_array[i].fyy,
--                                               neighbour_init.padded_particle_array[i].fzz)
--        end
--    end
    klock = klock + 1
    --dl_meso_timing_mod.timchk(config) --TODO Include this.
    --TODO: NYI    frctim = timelp - timfst
    --TODO: NYI    timfrc = timfrc + frctim
    if variables.config[0].btype == 0 then
        [mdvv_nvt(2)];
    end

    --TODO: statis with some parallelism.
	dl_meso_statistics_mod.statis(variables.config, neighbour_init.padded_particle_array)

    --dl_meso_timing_mod.timchk(config) --TODO Include this.

    --TODO NYI: write_output_summary
    if variables.config[0].nstep % variables.config[0].nsbpo == 0 then
	    write_output_summary(0.0, false, false, variables.config);
        if variables.config[0].nstep == variables.config[0].nseql then
	    	write_output_equil(variables.config);
        end
    end

    --dl_meso_timing_mod.timchk(config) --TODO Include this.
--        stptim = timelp - timsrt
--        timstp = timstp + stptim
    variables.config[0].time = variables.config[0].tstep * double(variables.config[0].nstep - variables.config[0].nseql)
    if variables.config[0].nstep >= variables.config[0].nrun then --.OR. (timelp>=timjob-tclose)
        format.println("Finishing after step {} of {}", variables.config[0].nstep, variables.config[0].nrun)
        finish = true
    end

    var time : double = variables.config[0].tstep * double(variables.config[0].nstep - variables.config[0].nseql)
    if variables.config[0].lcorr and (variables.config[0].nstep >= variables.config[0].nseql ) and ( (variables.config[0].nstep - variables.config[0].nseql) % variables.config[0].iscorr == 0) then
        write_correl(variables.config, time)
    end

    --TODO NYI: write_stress

    l_hist = (variables.config[0].ltraj and (variables.config[0].nstep>=variables.config[0].straj) and ( (variables.config[0].nstep - variables.config[0].straj) % variables.config[0].ntraj == 0))
    l_exp = (not finish) and ((variables.config[0].nstep % variables.config[0].ndump) == 0)
    --TODO NYI: Gather write data and write revive
  
--    if l_hist or l_exp then
--        gather_write_data(l_exp, l_hist, time, variables.config, neighbour_init.padded_particle_array)
--    end
--    if l_exp then
--        write_revive(variables.config)
--   end
end
--***********************************************************************
--     end of dissipative particle dynamics calculations
--***********************************************************************
    variables.config[0].timfrc = variables.config[0].timfrc / double(klock)
    variables.config[0].timstp = variables.config[0].timstp / double(klock)
    --TODO NYI
--      IF (idnode==0) WRITE (nprint, "(/,1x,'run terminating. elapsed cpu time = ',f10.2, ', job time = ',f10.2,&
--                                   &', close time = ',f10.2,/)") timelp, timjob, tclose
    var OUTPUT = c_stdio.fopen('OUTPUT', 'a')    
    c_stdio.fprintf(OUTPUT, "\n run terminating. elapsed cpu time = %10.2f, job time = %10.2f, close time = %10.2f\n",
                    0.0, 0.0, 0.0)
    --c_stdio.fprintf(OUTPUT, "\n run terminating. elapsed cpu time = %10.2f, job time = %10.2f, close time = %10.2f\n",
    --                variables.config[0].timelp, variables.config[0].timjob, variables.config[0].tclose)
    c_stdio.fclose(OUTPUT)
    write_output_result(variables.config, neighbour_init.padded_particle_array)
    --dl_meso_timing_mod.timchk(config) --TODO Include this.
end

--run_DSL(mdvv)
regentlib.start(mdvv)

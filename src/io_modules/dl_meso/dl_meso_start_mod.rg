import "regent"

local c_stdio = terralib.includec("stdio.h")

local c = regentlib.c

local fkt = 2.0/3.0

local start_mod = {}

--Finds the lowest value
local terra min_loc( array : int32[4]) : int
    var minval : int32 = 2147483647
    var minpos = -1
    var i = 0
    for i = 0, 4 do
        if array[i] < minval then
            minpos = i
            minval = array[i]
        end
    end
    return minpos
end

local cbrt = regentlib.cbrt(double)
local sqrt = regentlib.sqrt(double)

local task initialize(parts :  region(ispace(int1d), part), config : region(ispace(int1d), config_type)) where writes(parts), reads(parts), reads(config)
do

    --Initialise randomness --NYI: Use seed from CONTROL
    c.srand48(123)

    var gblnbd : &int32
    var gbltotnbd : &int32
    var totnbd : &int32

    --Allocate arrays for determining distribution of unbonded particles
--    gblnbd = [&&int32](regentlib.c.malloc([terralib.sizeof(&int32)] * 1))
    gblnbd = [&int32](regentlib.c.malloc([terralib.sizeof(int32)] * config[0].nspe))
    gbltotnbd = [&int32](regentlib.c.malloc([terralib.sizeof(int32)] * 1))
    totnbd = [&int32](regentlib.c.malloc([terralib.sizeof(int32)] * 1))

--construct lattices of beads and insert molecules

--insert frozen walls as simple cubic lattices (TODO: NYI lines 373-486)

--calculate numbers for orthorhombic lattice of particles (for all unbonded particles
--apart from those in frozen bead walls)

    var nplx : int = 0
    var nply : int = 0
    var nplz : int = 0
    var nplxshift : int = 0
    var nplyshift : int = 0
    var nplzshift : int = 0

    --TODO: This is nusyst - nfwsyst
    if config[0].nusyst > 0 then
        var yxratio : double
        var zxratio : double
        --TODO: These values actually depend on frzw*wid but we don't use this yet.
        yxratio = config[0].space.dim_y / config[0].space.dim_x
        zxratio = config[0].space.dim_z / config[0].space.dim_x

        --find number of beads required in each direction (TODO: also should be nusyst-nfwsyst)
        var nbds : double = cbrt(double(config[0].nusyst) / (yxratio * zxratio) )
        var nbdx : int32 = int32(nbds)
        var nbdy : int32 = int32(nbds*yxratio)
        var nbdsint : int32
        var xs : int32[4]
        xs[0] = (config[0].nusyst) % (nbdx*nbdy)
        if xs[0] > 0 then xs[0] = nbdx*nbdy - xs[0] end
        xs[1] = (config[0].nusyst) % ((nbdx+1)*nbdy)
        if xs[1] > 0 then xs[1] = (nbdx+1)*nbdy - xs[1] end
        xs[2] = (config[0].nusyst) % (nbdx*(nbdy+1))
        if xs[2] > 0 then xs[2] = nbdx*(nbdy+1) - xs[2] end
        xs[3] = (config[0].nusyst) % ((nbdx+1)*(nbdy+1))
        if xs[2] > 0 then xs[3] = (nbdx+1)*(nbdy+1) - xs[3] end
        var xsindex = min_loc(xs)
        
        if xsindex == 0 then
            nbdsint = xs[0]
        elseif xsindex == 1 then
            nbdsint = xs[1]
            nbdx = nbdx + 1
        elseif xsindex == 2 then
            nbdsint = xs[2]
            nbdy = nbdx + 1
        elseif xsindex == 3 then
            nbdsint = xs[3]
            nbdx = nbdx + 1
            nbdy = nbdy + 1
        end
        var nbdz : int32 = (config[0].nusyst + nbdsint) / (nbdx * nbdy)

        --calculate spacings for orthorhombic lattice (TODO also should be space.dim - 2.0*frzwid)
        var disx : double = (config[0].space.dim_x) / double(nbdx) 
        var disy : double = (config[0].space.dim_y) / double(nbdy) 
        var disz : double = (config[0].space.dim_z) / double(nbdz) 
      
        --Don't worry about calculating spaces for based on MPI
        --So set values accordingly
        nplx = nbdx
        nply = nbdy
        nplz = nbdz
        var xdispl : double = 0.5 * disx
        var ydispl : double = 0.5 * disy
        var zdispl : double = 0.5 * disz

        var nubeads : int1d = int1d(0)
        var shfz : double = zdispl
        for k = 0, nplz do
            var shfy = ydispl
            for j = 0, nply do
                var shfx = xdispl
                for i = 0, nplx do
                    var isx = i 
                    var isy = j
                    var isz = k
                    --TODO: should include nfsyst
                    var lb = 1 + isx  + nbdx*( isy + nbdy * isz)
                    
                    if lb <= config[0].nusyst then
                        parts[nubeads].core_part_space.pos_x = shfx
                        parts[nubeads].core_part_space.pos_y = shfy
                        parts[nubeads].core_part_space.pos_z = shfz
                        parts[nubeads].lab = lb
                        parts[nubeads].ltm = 0
                        nubeads = nubeads + int1d(1)
                    end
                    shfx = shfx + disx
                end
                shfy = shfy + disy
            end
            shfz = shfz + disz
        end
    end

    --Put in bonded beads (TODO NYI) lines 618-790

    --Check total number of beads TODO NYI lines 794-801

    --Assign species to unbonded particles
    var nuspe = 0
    var ispe = 0
    for i=0, config[0].nspe do
        if config[0].nspec[i] > 0 then
            nuspe = nuspe + 1
            ispe = i
        end
    end
    
    if nuspe == 1 then
        for j in parts.ispace do
            parts[j].ltp = ispe
        end
    elseif nuspe > 1 then
    --randomize unbonded particle assignments
        var gblnbd : &int32 = [&int32](regentlib.c.malloc([terralib.sizeof(int32)] * config[0].nspe))
        var totnbd : int32 = config[0].nusyst
        for ispec = 0, config[0].nspe do
            --TODO NYI remove frozen from count when implemented
            var specconc : double = (double(config[0].nspec[ispec]) / double(config[0].nusyst))
            gblnbd[ispec] = int32( double(config[0].nusyst) * specconc)
            totnbd = totnbd - gblnbd[ispec]
        end
        --TODO NYI: Frozen specie management

        ispe = -1
        var ntop = 0
        for i in parts.ispace do
            while int32(i) >= ntop do
                ispe = ispe + 1
                ntop = ntop + gblnbd[ispe]
            end
            parts[i].ltp = ispe
        end

        ntop = regentlib.fmax(config[0].nusyst/20, 100)
        ntop = 0
        --Performan random swaps
        for j = 0, ntop do
            --Using inbuilt randomness for now?
            for k = 0, config[0].nusyst do
                var temp : int32 = int32( double(config[0].nusyst) * c.drand48())

                var temp_ltp = parts[int1d(k)].ltp
                parts[int1d(k)].ltp = parts[int1d(temp)].ltp
                parts[int1d(temp)].ltp = temp_ltp 
            end
        end

    end

    --TODO NYI: if frozen beads call sort beads line 898

    regentlib.c.free(gbltotnbd)
    regentlib.c.free(totnbd)
--    regentlib.c.free(gblnbd[0])
    regentlib.c.free(gblnbd)

end

task start_mod.initialvelocity(parts : region(ispace(int1d), part), config :region(ispace(int1d), config_type)) where writes(parts), reads(parts), reads(config)
do
-- determine system total mass (omitting frozen beads)
    var syswgt : double = 0.0
    for i = 0, config[0].nspe do
        --Don't count frozen beads TODO: NYI
        syswgt = syswgt + config[0].masstmp[i] * double(config[0].nspec[i])
    end

    --Reset forces to zero
    for i in parts.ispace do
        parts[i].fxx = 0.0
        parts[i].fyy = 0.0
        parts[i].fzz = 0.0
        if config[0].lvarfc then
            parts[i].fvx = 0.0
            parts[i].fvy = 0.0
            parts[i].fvz = 0.0
        end
    end

    --Set starting velocities
    var buf : double[6]
    buf[0] = 0.0
    buf[1] = 0.0
    buf[2] = 0.0
    buf[3] = 0.0
    buf[4] = 0.0
    buf[5] = 0.0
    --TODO NYI: Freeze frozen beads
    for i in parts.ispace do
        parts[i].core_part_space.vel_x = c.drand48() - 0.5
        parts[i].core_part_space.vel_y = c.drand48() - 0.5
        parts[i].core_part_space.vel_z = c.drand48() - 0.5
        buf[0] = buf[0] + parts[i].core_part_space.mass * parts[i].core_part_space.vel_x
        buf[1] = buf[1] + parts[i].core_part_space.mass * parts[i].core_part_space.vel_y
        buf[2] = buf[2] + parts[i].core_part_space.mass * parts[i].core_part_space.vel_z
        buf[3] = buf[3] + parts[i].core_part_space.mass * parts[i].core_part_space.vel_x * parts[i].core_part_space.vel_x
        buf[4] = buf[4] + parts[i].core_part_space.mass * parts[i].core_part_space.vel_y * parts[i].core_part_space.vel_y
        buf[5] = buf[5] + parts[i].core_part_space.mass * parts[i].core_part_space.vel_z * parts[i].core_part_space.vel_z
    end

    buf[0] = buf[0] / syswgt
    buf[1] = buf[1] / syswgt
    buf[2] = buf[2] / syswgt
    buf[3] = 0.5 * fkt * (buf[3] - syswgt * buf[0] * buf[0]) / double(config[0].nsyst)
    buf[4] = 0.5 * fkt * (buf[4] - syswgt * buf[1] * buf[1]) / double(config[0].nsyst)
    buf[5] = 0.5 * fkt * (buf[5] - syswgt * buf[2] * buf[2]) / double(config[0].nsyst)
    var tscal : double = sqrt(config[0].temp / (buf[3] + buf[4] + buf[5]))
    --TODO :NYI Ignore frozen beads
    var maxx = -1.0
    var maxy = -1.0
    var maxz = -1.0
    for i in parts.ispace do
        parts[i].core_part_space.vel_x = tscal * (parts[i].core_part_space.vel_x - buf[0])
        parts[i].core_part_space.vel_y = tscal * (parts[i].core_part_space.vel_y - buf[1])
        parts[i].core_part_space.vel_z = tscal * (parts[i].core_part_space.vel_z - buf[2])
    end
end

task start_mod.start( parts : region(ispace(int1d), part), config :region(ispace(int1d), config_type)) where writes(parts), reads(parts), reads(config),
                                                                            writes(config)  --Overapproximating config writes to avoid Legion issue #289
do

--TODO NYI: set restart filename

--TODO: NYI allocate arrays for bonds

--TODO: NYI set up array for global bead numbers of molecules (for unit cell)

--TODO: NYI allocate global/local list for bonded particles

--restart control option
    if config[0].kres == 0 then
        if config[0].l_config then
            regentlib.assert(false, "CONFIG files not yet supported")
        else
            initialize(parts, config)
        end
    elseif config[0].kres == 1 then
        --restart control option: restart from end point of previous simulation
        regentlib.assert(false, "Restarting not yet supported")
    elseif config[0].kres == 2 then
        --restart control option: restart from end point as new simulation, no temperature reset
        regentlib.assert(false, "Restarting not yet supported")
    elseif config[0].kres == 3 then
        --restart control option: restart from end point as new simulation with temperature reset
        regentlib.assert(false, "Restarting not yet supported")
    end

--TODO: NYI nfold not supported lines 151-184

--TODO: NYI assign bonds, constraints, angles and dihedrals to tables

-- print the link cell numbers per domain FIXME: printing 1 for now...
    var OUTPUT = c_stdio.fopen("OUTPUT", "a")
    c_stdio.fprintf(OUTPUT, "\n link cell details\n      number of cells per domain (npx,npy,npz)=     1    1    1\n")
    --if etype > 0 then TODO: NYI print electrostatic link cell details
    --if lmb > 0 then TODO: NYI print manybodylink cell details

--Assign particle/molecule names and masses
    for i in parts.ispace do
        parts[i].core_part_space.mass = config[0].masstmp[parts[i].ltp]
        for j =0, 9 do
            parts[i].atmnam[j] = config[0].namspe[parts[i].ltp][j]
            --TODO: NYI parts[i].molnam[j] = config[0].mol[parts[i].ltm][j]
        end
    end

    --assign particle velocities to give required system temperature
    if config[0].kres == 0 then --TODO config NYIand ( (not config[0].l_config) or (config[0].l_config and config[0].levcfg == 0)) then
        start_mod.initialvelocity(parts, config)
    end

--TODO: NYI quench velocities for particles with constraints lines 218-241

--TODO: NYI allocate pair lists for velocity correction thermostats lines 245-251

--TODO: NYI lines 257-264
--   !     determine (initial) reciprocal image search pattern and
--   !     calculate corrections to forces, potential energy etc. between frozen
--   !     charged beads (only need to do these once for NVT ensembles) 

--write out initial positions and velocities
    if config[0].nsyst > 0 then
        c_stdio.fprintf(OUTPUT, "\n initial particle positions and velocities\n")
        var k = regentlib.fmax(1, config[0].nsyst/20)
        var until_next_print :int32 = 0
        for i in parts.ispace do
            if until_next_print == int32(0) then
                until_next_print = k-1
                --if bond/cons then
                --else
                    var buffer : int8[25]
                    c_stdio.snprintf(&(buffer[0]), 25, "%i", parts[i].lab)
                    c_stdio.fprintf(OUTPUT, " %10s    %8s    %14.6e    %14.6e    %14.6e    %14.6e    %14.6e    %14.6e\n", buffer, parts[i].atmnam,
                                            parts[i].core_part_space.pos_x, parts[i].core_part_space.pos_y, parts[i].core_part_space.pos_z,
                                            parts[i].core_part_space.vel_x, parts[i].core_part_space.vel_y, parts[i].core_part_space.vel_z)
                --end
            else
                until_next_print = until_next_print - 1
            end
        end
    end
    c_stdio.fclose(OUTPUT)

    if config[0].ltraj then
        --Ignore wtype TODO : NYI
        [dl_meso_write_mod.write_history_header](config, parts)
    end
end

return start_mod

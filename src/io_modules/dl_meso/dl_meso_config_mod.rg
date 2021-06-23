import "regent"

local config_mod = {}

local c_unistd = terralib.includec("unistd.h")
local c_stdio = terralib.includec("stdio.h")

local sqrt = regentlib.sqrt(double)

task config_mod.sysdef(config : region(ispace(int1d), config_type)) where writes(config), reads(config)
do

--    config[0].nummol = 0
--    config[0].nummolcell = 0
--  config[0].nbonddef = 0
--  config[0].ncondef = 0
--  config[0].nangdef = 0
--  config[0].ndhddef = 0
--  config[0].nfwsyst = 0
--    config[0].itype = 0
--    config[0].btype = 0
--  config[0].etype = 0
--  config[0].srftype = 0
--  config[0].srfx = 0
--  config[0].srfy = 0
--  config[0].srfz = 0
--  config[0].nfoldx = 1
--  config[0].nfoldy = 1
--  config[0].nfoldz = 1
--  config[0].nlx = 0
--  config[0].nly = 0
--  config[0].nlz = 0
--  config[0].nlx2 = 0
--  config[0].nly2 = 0
--  config[0].nlz2 = 0
--  config[0].nlewx = 0
--  config[0].nlewy = 0
--  config[0].nlewz = 0
--  config[0].nlewx2 = 0
--  config[0].nlewy2 = 0
--  config[0].nlewz2 = 0
--  config[0].nlmbx = 0
--  config[0].nlmby = 0
--  config[0].nlmbz = 0
--  config[0].nlmbx2 = 0
--  config[0].nlmby2 = 0
--  config[0].nlmbz2 = 0
--  config[0].kmax1 = 0
--  config[0].kmax2 = 0
--  config[0].kmax3 = 0
  config[0].ndump = 1000
--  config[0].bdfrcx = 0.0
--  config[0].bdfrcy = 0.0
--  config[0].bdfrcz = 0.0
--  config[0].shrvx = 0.0
--  config[0].shrvy = 0.0
--  config[0].shrvz = 0.0
--  config[0].elecx = 0.0
--  config[0].elecy = 0.0
--  config[0].elecz = 0.0
--  config[0].vgapx = 0.0
--  config[0].vgapy = 0.0
--  config[0].vgapz = 0.0
--  config[0].dvar = 0.0
--    config[0].cutoff = 0.0 --rcut
--  config[0].rtcut = 0.0
--  config[0].rhalo = 0.0
--  config[0].relec = 0.0
--  config[0].rmbcut = 0.0
--  config[0].srfzcut = 0.0
--  config[0].srfpos = 0.0
 --   config[0].temp = 0.0
 --   config[0].tstep = 0.0
--    config[0].rtstep = 0.0
--    config[0].space.dim_x = 0.0 --dimx
--    config[0].space.dim_y = 0.0 --dimy
--    config[0].space.dim_z = 0.0 --dimz
--  config[0].dimxcell = 0.0
--  config[0].dimycell = 0.0
--  config[0].dimzcell = 0.0
--  config[0].gammaelec = 0.0
--  config[0].bjerelec = 0.0
--  config[0].alphaew = 0.0
--  config[0].betaew = 0.0
--  config[0].chglen = 0.0
--  config[0].delpot = 0.0
--  config[0].cutpot = 0.0
--  config[0].potrdr = 0.0
--  config[0].mxitshake = 250
--  config[0].shaketol = 1e-6
--    config[0].nstep = 0
--    config[0].nseql = 0
--  config[0].nshrs = 0
--    config[0].straj = 0
--    config[0].ntraj = 0
--    config[0].keytrj = 0
--  config[0].sstrs = 0
--  config[0].nstrs = 0
--  config[0].engunit = 0
--  config[0].rndseed = 0
--  config[0].mxspl = 0
--  config[0].npotgrid = 0
--  config[0].gridnpot = 0
--  config[0].wtype = 0
--  config[0].nminwriter = 4
--  config[0].gathermem = 32.0
--  config[0].lbond = false
--  config[0].lcons = false
--  config[0].langle = false
--  config[0].ldihed = false
--  config[0].lgbnd = false
--  config[0].lmb = false
--  config[0].ldpol = false
--  config[0].lnfold = false
--    config[0].ldyn = false
--  config[0].ligindex = false
--    config[0].ltraj = false
--  config[0].lstrs = false
--  config[0].lstrs[0:3] = false -- Split
--  config[0].lfrzwall = false
--  config[0].lfrzx = false
--  config[0].lfrzy = false
--  config[0].lfrzz = false
--  config[0].lisoprs = false
--   config[0].lcorr = false
--    config[0].ltemp = false
--  config[0].lvarfc = false
--  config[0].lconfzero = false
--  config[0].lompcrit = false
--    config[0].ltabpot = false

--TODO: NYI Config files not yet supported
    if( (c_unistd.access("CONFIG", c_unistd.R_OK) == 0 ) and config[0].l_config) then
        config[0].l_config = true
        regentlib.assert(false, "Reading CONFIG files not yet supported")
        --scan_config(...)
    end

    --TODO: NYI Check for restarting information
    if config[0].l_rest then
        --scan_export(...)
        regentlib.assert(false, "Reading export files not yet supported")
    end


--  read in system parameters and job title
    --First argument is : (not ((l_config and imcon > 0) or l_rest))
    --Since l_config and l_rest both must be false then this is true for now.
    dl_meso_read_mod.read_control(config, true, config[0].l_config, config[0].l_rest)

    --Write the name to the OUTPUT file
    var OUTPUT = c_stdio.fopen('OUTPUT', 'a')
    c_stdio.fprintf(OUTPUT,"%s\n", [rawstring](config[0].text))
    c_stdio.fclose(OUTPUT)

    dl_meso_read_mod.scan_field(config)

    --TODO: NYI table support
    --if(ltabpot) scan_table(...)

    dl_meso_read_mod.read_field(config)

    if config[0].cutoff < 1e-16 then
        config[0].cutoff = config[0].rtcut
    end 
    regentlib.assert(config[0].cutoff >= 1e-16, "Error 1 in read_control, rcut too small")
    if config[0].rtcut < 1e-16 then
        config[0].rtcut = config[0].cutoff
    end 
    config[0].rtct2 = config[0].rtcut * config[0].rtcut
    config[0].rrtcut = 1.0 / config[0].rtcut

--TODO: NYI - no support for these yet
    --if etype ... then relec = rcut end
--  config[0].rel2 = relec * relec
    --if config[0].srftype > 1 then
    --end
--srfzct2 = srfzcut * srfzcut

--      IF (lmb .AND. (rmbcut<1.0e-16_dp .OR. rmbcut>rcut)) rmbcut = rcut
--      rmbct2 = rmbcut * rmbcut
--      rrmbcut = MERGE (1.0_dp / rmbcut, 0.0_dp, lmb)
--
--      rcut = MAX (rcut, rtcut, srfzcut)

    config[0].rct2 = config[0].cutoff * config[0].cutoff
    config[0].rrct2 = 1.0 / config[0].rct2

    --TODO NYI: re-assign halo size according to maximum possible bond or constraint length (line 205-229)

    --TODO NYI: set tolerances for constraints (line 232)

    -- determine additional thermostat parameters
    if config[0].itype >= 0 and config[0].itype <= 3 then
        for i = 0, config[0].npot do
            config[0].sigma[i] = sqrt(2.0 * config[0].gamma[i] * config[0].temp * config[0].rtstep)
        end
    elseif config[0].itype == 4 then
        for i = 0, config[0].npot do
            config[0].sigma[i] = config[0].gamma[i] * config[0].tstep
        end
    elseif config[0].itype == 5 then
        --Nothing here in dl_meso
    elseif config[0].itype == 6 then
        for i = 0, config[0].npot do
            config[0].sigma[i] = config[0].gamma[i] * config[0].tstep
        end
    end

    if ((config[0].itype > 0 and config[0].itype < 4) or config[0].itype == 6) then
        config[0].lvarfc = true
    else
        config[0].lvarfc = false
    end

    if config[0].btype > 0 and config[0].btype <= 3 then
        config[0].psmass = double(config[0].nsyst-config[0].nfsyst) * config[0].temp * config[0].abaro * config[0].abaro
        config[0].rpsmass = 1.0 / config[0].psmass
        config[0].sigmalang = sqrt(config[0].fkt * config[0].bbaro * config[0].psmass * config[0].temp * config[0].rtstep)
    end

    --TODO: determines number of additional frozen beads for walls and change to system volume (line 266 to 282)
    --TODO: determines position of hard surface or if shear can actually be applied (line 286-291)
    --select statistical properties to print
    config[0].outsel = 0
    --TODO: Add support for other outputs (lines 296-303)

    --No domain decomposition here.
    --TODO: NYI Support determine nodes with walls/surfaces

    --Calculate domain dimensions - think unneeded.
    --TODO NYI elecgen not yet available
    --TODO: NYI determine parameters for electrostatics

    --Ignore halo sizes for now

    --Ignore maximum array extents for now (regions/partitions should handle)
    --No explicit transfer buffers
    --TODO: NYI determine how many export particle lists to store - may not need
    --TODO: NYI particle pairs for non-DPD thermostats - may not need

    -- Don't think i need: check number of beads in cell domain
    -- TODO: Don't think needed -  Setup output file groups

    OUTPUT = c_stdio.fopen('OUTPUT', 'a')
    --write out system data
    c_stdio.fprintf(OUTPUT, "\n physical specification\n\n")
    c_stdio.fprintf(OUTPUT, "      system volume                 = %14.6e\n", config[0].volm)
    c_stdio.fprintf(OUTPUT, "      simulation cell dimensions    = %14.6e  %14.6e  %14.6e\n", config[0].space.dim_x,
                                                                config[0].space.dim_y, config[0].space.dim_z)
    --if config[0].nfold then
    --end
    c_stdio.fprintf(OUTPUT, "      specified temperature         = %14.6e\n", config[0].temp)
    --if config[0].btype > 0 then
    --    c_stdio.fprintf(OUTPUT, "      specified pressure            = %14.6e\n", config[0].prszero)
    --end
    c_stdio.fprintf(OUTPUT, "      simulation timestep           = %14.6e\n", config[0].tstep)

    c_stdio.fprintf(OUTPUT, "\n simulation controls\n\n")
    if config[0].kres == 0 then
        c_stdio.fprintf(OUTPUT, "      new simulation\n")
    elseif config[0].kres == 1 then
        c_stdio.fprintf(OUTPUT, "      restarting simulation\n")
    elseif config[0].kres == 2 then
        c_stdio.fprintf(OUTPUT, "      new simulation with initial state from restart files\n")
    elseif config[0].kres == 3 then
        c_stdio.fprintf(OUTPUT, "      new simulation with initial state from restart files and temperature reset\n")
    end

    var buffer : int8[25]
    c_stdio.snprintf(&(buffer[0]), 25, "%i", config[0].nsyst)
    c_stdio.fprintf(OUTPUT, "      system size (particles)      = %14s\n", buffer)
    
    c_stdio.snprintf(&(buffer[0]), 25, "%i", 0) --TODO: num frozen parts
    c_stdio.fprintf(OUTPUT, "      number of frozen particles   = %14s\n", buffer)

    c_stdio.snprintf(&(buffer[0]), 25, "%i", config[0].nsyst)
    c_stdio.fprintf(OUTPUT, "      max no. of particles/node    = %14s\n", buffer)

    c_stdio.snprintf(&(buffer[0]), 25, "%i", 0) --TODO: buffer parts
    c_stdio.fprintf(OUTPUT, "      init. no. buffer particles   = %14s\n", buffer)
    
    c_stdio.snprintf(&(buffer[0]), 25, "%i", config[0].nrun)
    c_stdio.fprintf(OUTPUT, "      number of timesteps          = %14s\n", buffer)

    c_stdio.snprintf(&(buffer[0]), 25, "%i", config[0].ndump)
    c_stdio.fprintf(OUTPUT, "      restart file interval        = %14s\n", buffer)

    c_stdio.snprintf(&(buffer[0]), 25, "%i", config[0].nsbpo)
    c_stdio.fprintf(OUTPUT, "      printing interval            = %14s\n", buffer)

    c_stdio.snprintf(&(buffer[0]), 25, "%i", config[0].ltraj)
    c_stdio.fprintf(OUTPUT, "      data saving option           = %14s\n", buffer)

    c_stdio.snprintf(&(buffer[0]), 25, "%i", config[0].straj)
    c_stdio.fprintf(OUTPUT, "      data saving starting step    = %14s\n", buffer)

    c_stdio.snprintf(&(buffer[0]), 25, "%i", config[0].ntraj)
    c_stdio.fprintf(OUTPUT, "      data saving interval         = %14s\n", buffer)

    c_stdio.snprintf(&(buffer[0]), 25, "%i", config[0].keytrj)
    c_stdio.fprintf(OUTPUT, "      data saving information key  = %14s\n", buffer)

    c_stdio.snprintf(&(buffer[0]), 25, "%i", config[0].nstk)
    c_stdio.fprintf(OUTPUT, "      data stacking interval       = %14s\n", buffer)

    if config[0].ltemp then
        c_stdio.snprintf(&(buffer[0]), 25, "T")
    else
        c_stdio.snprintf(&(buffer[0]), 25, "F")
    end
    c_stdio.fprintf(OUTPUT, "      temperature scaling option   = %14s\n", buffer)

    c_stdio.snprintf(&(buffer[0]), 25, "%i", config[0].nsbts)
    c_stdio.fprintf(OUTPUT, "      temperature scaling interval = %14s\n", buffer)

    c_stdio.snprintf(&(buffer[0]), 25, "%i", config[0].nseql)
    c_stdio.fprintf(OUTPUT, "      equilibration period         = %14s\n", buffer)


--  write information about trajectory file format and writing properties - TODO: these are just 1 for now
    c_stdio.snprintf(&(buffer[0]), 25, "%i", 1)
    c_stdio.fprintf(OUTPUT, "\n      file writing options\n\n")

    c_stdio.fprintf(OUTPUT, "      min. no. of file writers     = %14s\n", buffer)
    c_stdio.fprintf(OUTPUT, "      actual no. of file writers   = %14s\n", buffer)
    c_stdio.fprintf(OUTPUT, "      max. gathered data (MiB)     = %14.6f\n", 0.0)

    --What is data saving format (Only DL_MESO supported)
    c_stdio.fprintf(OUTPUT, "\n      data saving format: standard (native DL_MESO_DPD)\n")

    --write selection and parameters for separated tensor files - TODO: NYI
    --if lstrs then
    --end

--!     write out cutoff and halo sizes (warn if latter are larger than half size of subdomain,
--!     abort if using SPME and halo sizes are larger than subdomain)
    c_stdio.fprintf(OUTPUT, "\n      maximum cutoff radius        = %14.6e\n", config[0].cutoff)
    c_stdio.fprintf(OUTPUT, "\n      thermostat cutoff radius     = %14.6e\n", config[0].rtcut)
    --TODO Not yet supported:
    --if lmb ...
    --if srftype > 0 ...
    --if etype > 0 ...
    c_stdio.fprintf(OUTPUT, "      domain boundary halo size    = %14.6e\n", 0.0) -- TODO: NYI slash ignored

    --if rhalox > 0.5*sidex .... abort
    -- write out integrator and barostat types and parameters
    if config[0].itype == 0 then
        c_stdio.fprintf(OUTPUT, "\n      integrator/thermostat type: md velocity verlet\n")
    elseif config[0].itype == 1 then
        c_stdio.fprintf(OUTPUT, "\n      integrator/thermostat type: dpd velocity verlet\n")
    elseif config[0].itype == 2 then
        c_stdio.fprintf(OUTPUT, "\n      integrator/thermostat type: dpd with first-order shardlow splitting\n")
    elseif config[0].itype == 3 then
        c_stdio.fprintf(OUTPUT, "\n      integrator/thermostat type: dpd with second-order shardlow splitting\n")
    elseif config[0].itype == 4 then
        c_stdio.fprintf(OUTPUT, "\n      integrator/thermostat type: lowe-andersen\n")
    elseif config[0].itype == 5 then
        c_stdio.fprintf(OUTPUT, "\n      integrator/thermostat type: peters\n")
    elseif config[0].itype == 6 then
        c_stdio.fprintf(OUTPUT, "\n      integrator/thermostat type: stoyanov-groot (lowe-andersen/nose-hoover)\n")
        --c_stdio.fprintf(OUTPUT, "      nose-hoover coupling param.  = %14.6e\n", config[0].alphasg)
        regentlib.assert(false, "Not yet supporting integrator type stoyanov-groot")
    end

    if config[0].btype == 0 then
        c_stdio.fprintf(OUTPUT, "\n      ensemble type: NVT (constant volume/temperature)\n")
        c_stdio.fprintf(OUTPUT, "      barostat type: none\n")
    elseif config[0].btype == 1 then
        regentlib.assert(false, "Only currently supporting NVT (constant volume/temperature) with no barostat")
    elseif config[0].btype == 2 then
        regentlib.assert(false, "Only currently supporting NVT (constant volume/temperature) with no barostat")
    elseif config[0].btype == 3 then
        regentlib.assert(false, "Only currently supporting NVT (constant volume/temperature) with no barostat")
    elseif config[0].btype == 4 then
        regentlib.assert(false, "Only currently supporting NVT (constant volume/temperature) with no barostat")
    elseif config[0].btype == 5 then
        regentlib.assert(false, "Only currently supporting NVT (constant volume/temperature) with no barostat")
    elseif config[0].btype == 6 then
        regentlib.assert(false, "Only currently supporting NVT (constant volume/temperature) with no barostat")
    end

    c_stdio.fprintf(OUTPUT, "\n species data\n")
    c_stdio.fprintf(OUTPUT, "\n                population        mass            charge            frozen\n")

    var qisyst : double = 0.0
    for i=0, config[0].nspe do
        --TODO: Nspec not yet implemented FIXME: URGENT READ_FIELD
        c_stdio.snprintf(&(buffer[0]), 25, "%i", config[0].nspec[i])
        c_stdio.fprintf(OUTPUT, "      %8s  %14s  %14.6e  %14.6e           F\n", config[0].namspe[i], buffer, config[0].masstmp[i], config[0].chgetmp[i])
        qisyst = qisyst + double(config[0].nspec[i]) * config[0].chgetmp[i] --TODO: nspec + nspecmol once molecules exist
        --NB last F is a %10s storing lfrzn(i) > 0 -- no frozen so false at the moment
    end
    if qisyst > 1e-10 or qisyst < -1e-10 then
        regentlib.assert(false, "Error -6 in this function....")
    end

    --TODO NYI: write out energy scaling lines 656-663

    --TODO NYI: write out randomisation seed (if supplied) lines 665-668

--  write out interaction data
    c_stdio.fprintf(OUTPUT, "\n energy parameters: absolute values\n")
    c_stdio.fprintf(OUTPUT, "\n interaction potential parameters\n")
    if config[0].ltabpot then
        regentlib.assert(false, "ltabpot not yet supported")
--        c_stdio.fprintf(OUTPUT, "\n      potential tables read from TABLE file\n")
--        c_stdio.fprintf(OUTPUT, "\n      table grid spacing           = %14.6e\n", config[0].delpot)
--        c_stdio.fprintf(OUTPUT, "\n      potential cutoff radius      = %14.6e\n", config[0].cutpot)
--        c_stdio.fprintf(OUTPUT, "\n      no. of grid points/potential = %14.6e\n", config[0].npotgrid)
    end
    
    c_stdio.fprintf(OUTPUT, "\n      energy parameters\n\n")
    for i=0, config[0].nspe do
        for j=i, config[0].nspe do
            var k = ((j+1) * j) / 2 + (i+1) - 1
            if config[0].ktype[k] == -1 then
                regentlib.assert(false, "ktype -1 not yet supported")
--                c_stdio.fprintf(OUTPUT, "      %8s %8s  tab  %14.6  %14.6e\n", config[0].namspe[i], config[0].namspe[j],
--                        MINVAL (tab_potential (1:(npotgrid+1), k)), MAXVAL (tab_potential (1:(npotgrid+1), k))
            elseif config[0].ktype[k] == 0 then
                c_stdio.fprintf(OUTPUT, "      %8s  %8s  lj    %14.6e\n", config[0].namspe[i], config[0].namspe[j], 0.25 * config[0].vvv[0][k])
            elseif config[0].ktype[k] == 1 then
                c_stdio.fprintf(OUTPUT, "      %8s  %8s  wca   %14.6e\n", config[0].namspe[i], config[0].namspe[j], 0.25 * config[0].vvv[0][k])
            elseif config[0].ktype[k] == 2 then
                c_stdio.fprintf(OUTPUT, "      %8s  %8s  dpd   %14.6e\n", config[0].namspe[i], config[0].namspe[j], config[0].vvv[0][k])
            elseif config[0].ktype[k] == 3 then
                c_stdio.fprintf(OUTPUT, "      %8s  %8s  mdpd  %14.6e  %14.6  %14.6e  %14.6  %14.6e\n", config[0].namspe[i], 
                                            config[0].namspe[j], config[0].vvv[0][k], config[0].vvv[1][k], config[0].vvv[2][k],
                                            config[0].vvv[3][k], config[0].vvv[4][k])
            end
        end
    end
--TODO: Not yet supporting lj or tabulated
--IF (ANY (ktype==-1) .OR. ANY(ktype==0)) WRITE (nprint, "(/,1x,'potential lr correction = ',1pe14.6,/,&
--                                 &1x,'virial lr correction    = ',1pe14.6)") clr (1)/volm, clr (2)/volm
    c_stdio.fprintf(OUTPUT, "\n      interaction lengths\n\n")
    var cutcheck : bool = false
    for i=0, config[0].nspe do
        for j=i, config[0].nspe do
            var k = ((j+1) * j) / 2 + (i+1) - 1
            if config[0].ktype[k] == -1 then
                c_stdio.fprintf(OUTPUT, "      %8s  %8s  %14.6e\n", config[0].namspe[i], config[0].namspe[j], config[0].vvv[0][k] )
                if config[0].vvv[0][k] > config[0].cutoff then
                    cutcheck = true
                end
            elseif config[0].ktype[k] == 0 then
                c_stdio.fprintf(OUTPUT, "      %8s  %8s  %14.6e\n", config[0].namspe[i], config[0].namspe[j], config[0].vvv[1][k] )
                if config[0].vvv[1][k] > config[0].cutoff then
                    cutcheck = true
                end
            elseif config[0].ktype[k] == 1 then
                c_stdio.fprintf(OUTPUT, "      %8s  %8s  %14.6e\n", config[0].namspe[i], config[0].namspe[j], config[0].vvv[1][k] )
                if config[0].vvv[1][k] > config[0].cutoff then
                    cutcheck = true
                end
            elseif config[0].ktype[k] == 2 then
                c_stdio.fprintf(OUTPUT, "      %8s  %8s  %14.6e\n", config[0].namspe[i], config[0].namspe[j], config[0].vvv[1][k] )
                if config[0].vvv[1][k] > config[0].cutoff then
                    cutcheck = true
                end
            elseif config[0].ktype[k] == 3 then
                c_stdio.fprintf(OUTPUT, "      %8s  %8s  %14.6e\n", config[0].namspe[i], config[0].namspe[j], config[0].vvv[5][k] )
                if config[0].vvv[5][k] > config[0].cutoff then
                    cutcheck = true
                end
            end
        end
    end
   
    regentlib.assert(cutcheck == false, "Error -7 in sysdef")

    --continue from 727 
    c_stdio.fprintf(OUTPUT, "\n      dissipative parameters\n")
    if config[0].itype >= 0 and config[0].itype <= 3 then
        c_stdio.fprintf(OUTPUT, "\n                                  viscosity       random force\n")
    elseif config[0].itype == 4 and config[0].itype == 6 then
        c_stdio.fprintf(OUTPUT, "\n                                  coll. freq.\n")
    elseif config[0].itype == 5 then
        c_stdio.fprintf(OUTPUT, "\n                                  viscosity\n")
    end
    for i=0, config[0].nspe do
        for j=i, config[0].nspe do
            var k = ((j+1) * j) / 2 + (i+1) - 1
            if config[0].itype >= 0 and config[0].itype <= 3 then
                c_stdio.fprintf(OUTPUT, "      %8s  %8s        %14.6e  %14.6e\n", config[0].namspe[i],  
                                            config[0].namspe[j], config[0].gamma[k], config[0].sigma[k])
            elseif config[0].itype >=4 and config[0].itype <= 6 then
                c_stdio.fprintf(OUTPUT, "      %8s  %8s        %14.6e\n", config[0].namspe[i],  
                                            config[0].namspe[j], config[0].gamma[k])
            end
        end
    end
   
    --TODO: NYI Frozen 
    -- IF (lfrzwall) THEN ... END (lines 748-758)
    --TODO: NYI srftype
    --SELECT CASE(srftype) ... (lines 760-819)

    --IF (ABS(bdfrcx >0.0... (lines 822-825)
    

--  TODO: NYI write out bond and molecule parameters (lines 827 to 911)

-- TODO: NYI   write out electrostatic parameters (lines 913-1063)

    regentlib.assert(config[0].l_scr == false, "Not yet supporting l_scr = true")
    c_stdio.fclose(OUTPUT)
    
--!     ensure history, correlation and stress tensor intervals are at
--!     least one to avoid floating-point exceptions during run
    config[0].ntraj = regentlib.fmax(config[0].ntraj, 1)
    config[0].iscorr = regentlib.fmax(config[0].iscorr, 1)
--    config[0].nstrs = regentlib.fmax(config[0].nstrs, 1)
--! TODO: NYI     prepare angle and dihedral parameters (convert to radians etc.) lines 1081-1105

--allocate arrays for bead positions, velocities, forces, potential energy etc. - TODO: Not done here - will be particle region

-- allocate transfer buffers only if running in parallel (bufnum=0 for serial version of code) or using SPME - Not needed i believe

--Not needed:
--!     allocate arrays to store particle indices during export communications
--!     (used to update velocities/positions for dpd-vv, shardlow or constraints)

--Not needed -already in config
--!     allocate statistical stack arrays for rolling averages

--Not needed !     assign pointers

--!Not needed:      if using OpenMP, check if using critical force assignments
--!     (based on user selection or available memory)
end


task config_mod.zero(config : region(ispace(int1d), config_type), parts : region(ispace(int1d), part)) where writes(config), reads(config), writes(parts), reads(parts)
do

    var x : double = 0.0
    
    --zero step counters
    config[0].nstep = 0
    config[0].nav = 0

    --initial time parameters
    config[0].timfrc = 0.0
    config[0].timstp = 0.0

    --set system parameters

    --config[0].pe = 0.0 --TODO: NYI
    --config[0].vir = 0.0 --TODO: NYI
    --config[0].tke = 0.0 --TODO: NYI
    --config[0].ee = 0.0 --TODO: NYI
    --config[0].be = 0.0 --TODO: NYI
    --config[0].ae = 0.0 --TODO: NYI
    --config[0].de = 0.0 --TODO: NYI
    --config[0].bdlng = 0.0 --TODO: NYO
    --config[0].bdlmin = 0.0 --TODO: NYI
    --config[0].bdlmax = 0.0 --TODO: NYI
    --config[0].bdang = 0.0 --TODO: NYI
    --config[0].bddhd = 0.0 --TODO: NYI
--    for i =0, 36 do
--        config[0].stress[i] = 0.0
--    end

    --set barostat parameters

    config[0].upx = 0.0
    config[0].upy = 0.0
    config[0].upz = 0.0
    config[0].fpx = 0.0
    config[0].fpy = 0.0
    config[0].fpz = 0.0

    --set shear parameters
    config[0].shrdx = 0.0
    config[0].shrdy = 0.0
    config[0].shrdz = 0.0

    --zero accumulators and stacks
    for i = 0, CONST_statsize do
        config[0].ave[i] = 0.0
        config[0].flc[i] = 0.0
    end
    for i = 0, 12 do
        config[0].zum[i] = 0.0
    end
    for i = 0, CONST_maxstk do
        config[0].stkpe[i] = 0.0  
        config[0].stkee[i] = 0.0  
        config[0].stkse[i] = 0.0  
        config[0].stkde[i] = 0.0  
        config[0].stkae[i] = 0.0  
        config[0].stkbe[i] = 0.0  
        config[0].stkvir[i] = 0.0 
        config[0].stkvlm[i] = 0.0
        config[0].stkzts[i] = 0.0 
        config[0].stktkex[i] = 0.0 
        config[0].stktkey[i] = 0.0 
        config[0].stktkez[i] = 0.0 
    end

    --zero velocities and forces
    for i in parts.ispace do
        parts[i].core_part_space.vel_x = 0.0
        parts[i].core_part_space.vel_y = 0.0
        parts[i].core_part_space.vel_z = 0.0
        parts[i].fxx = 0.0
        parts[i].fyy = 0.0
        parts[i].fzz = 0.0
        if config[0].lvarfc then
            parts[i].fvx = 0.0
            parts[i].fvy = 0.0
            parts[i].fvz = 0.0
        end
    end

    --TODO: NYI zero corrective values (electrostatics)


    --TODO: NYI start off random number generators...
end


return config_mod

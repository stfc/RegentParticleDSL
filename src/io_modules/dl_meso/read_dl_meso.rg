-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

local format = require("std/format")
local c_unistd = terralib.includec("unistd.h")
local c_stdlib = terralib.includec("stdlib.h")
local c_stdio = terralib.includec("stdio.h")
local c_string = terralib.includec("string.h")
local io_utils = require("src/io_modules/dl_meso/io_utils")
local c_math = terralib.includec("math.h")

local dl_meso_read_mod = {}

--Hack to deal with the fact that the EOF defition in stdio.h is not imported
--TODO: Improve this code because its linux-dependent right now...
--https://github.com/stfc/RegentParticleDSL/issues/39
local EOF = -1

task dl_meso_read_mod.scan_field --( cutoff : &double, srfzcut : double, mxprm : &int, ltabpot : bool, lrcut : bool, nspe : &int,
                --  namspe : &&&int8, masstmp : &&double, chgetmp : &&double )--,  ...)
               ( config : region(ispace(int1d), config_type)) where writes(config), reads(config) do
    var record : int8[201]
    var record1 : int8[201]

    var ifinish = 0

    var ipot : int
--    var imxprm : int = @mxprm
--    var rcut : double = @cutoff
    var inspe : int
 
    var key : &int8
    var word : &int8 = [&int8](regentlib.c.malloc(1))
  --Check FIELD file exists
  if( c_unistd.access("FIELD", c_unistd.R_OK) ~= 0 ) then
    format.println("FIELD file not found")
    c_stdlib.exit(1) 
  end

  --Open the FIELD file for first pass. From DL_MESO_DPD code:
  --determine numbers of species and
  --molecule types, maximum numbers of beads, bonds, constraints, angles
  --and dihedrals, maximum number of parameters needed for interaction
  --potentials and units for energy parameters
  var FIELD = c_stdio.fopen("FIELD", "r")
  --Skip the first line (this is a title)
  c_stdio.fscanf(FIELD,  "%*[^\n]\n")

  var finish : bool = false
  while (not finish) do
    --Read line
    var err = c_stdio.fgets(record, 200, FIELD)
    if err == [&int8](0) then
        finish = true
    end
    --Pull the first word
    key = io_utils.get_word(record, 1)
    --Convert to lowercase
    io_utils.to_lowercase(&key)

    if c_string.strncmp(key, "close", 5) == 0 then 
        finish = true
    end

    if c_string.strncmp(key, "finish", 6) == 0 then
        ifinish = ifinish + 1
    end

    if c_string.strncmp(key, "speci", 5) == 0 then
        regentlib.c.free(word)
        word = io_utils.get_word(record, 2)
        inspe = c_stdlib.atoi(word)
        config[0].nspe = inspe
        regentlib.assert( config[0].nspe < CONST_max_species, "Too many species in FIELD file, increase max_species")

    elseif c_string.strncmp(key, "molecul", 7) == 0 then
        --TODO: NYI
        regentlib.assert(false, "molecul keyword in FIELD not yet supported")
    
    elseif c_string.strncmp(key, "interact", 8) == 0 then
        --TODO: NYI
        regentlib.c.free(word)
        word = io_utils.get_word(record, 2)
        
        ipot = c_stdlib.atoi(word)

        for i = 1, ipot do
            c_stdio.fgets(record1, 200, FIELD)
            regentlib.c.free(word)
            word = io_utils.get_word(record1, 3)
            io_utils.to_lowercase(&word)
            if ( c_string.strncmp(word, "lj", 2) ) == 0 then
                config[0].mxprm = regentlib.fmax(config[0].mxprm, 3)
                if (config[0].lrcut) then config[0].cutoff = regentlib.fmax(config[0].cutoff, io_utils.get_double(record1, 5)) end
            elseif( c_string.strncmp(word, "wca", 3) ) == 0 then
                config[0].mxprm = regentlib.fmax(config[0].mxprm, 5)
                if (config[0].lrcut) then config[0].cutoff = regentlib.fmax(config[0].cutoff, io_utils.get_double(record1, 5)) end
            elseif( c_string.strncmp(word, "dpd", 3) ) == 0 then
                config[0].mxprm = regentlib.fmax(config[0].mxprm, 3)
                if (config[0].lrcut) then config[0].cutoff = regentlib.fmax(config[0].cutoff, io_utils.get_double(record1, 5)) end
            elseif (c_string.strncmp(word, "mdpd", 4) ) == 0 then 
                config[0].mxprm = regentlib.fmax(config[0].mxprm, 7)
                if (config[0].lrcut) then config[0].cutoff = regentlib.fmax(config[0].cutoff, io_utils.get_double(record1, 5)) end
            elseif( c_string.strncmp(word, "tab", 3) ) == 0 then
                config[0].mxprm = regentlib.fmax(config[0].mxprm, 2)
                config[0].ltabpot = true
            end
        end 


    elseif c_string.strncmp(key, "surf", 4 ) == 0 then
        --TODO: NYI


    elseif c_string.strncmp(key, "units", 5 ) == 0 then
        --TODO: NYI
    end 

    regentlib.c.free(key)

  end

  c_stdio.fclose(FIELD)

  config[0].npot = config[0].nspe * (config[0].nspe + 1) / 2
  --Allocate arrays
--  var temp_namspe : &&int8 = [&&int8](regentlib.c.malloc([terralib.sizeof(&int8)] * inspe))
--  for i = 0, inspe do
--    temp_namspe[i] = [&int8] (regentlib.c.malloc([terralib.sizeof(int8)]*9))
--  end
--  var temp_masstmp : &double = [&double] (regentlib.c.malloc([terralib.sizeof(double)] * inspe))
--  var temp_chgetmp : &double = [&double] (regentlib.c.malloc([terralib.sizeof(double)] * inspe))
   --Ignore bonds etc. for now
  --Second pass to pull species names (and in the future other things)
  FIELD = c_stdio.fopen("FIELD", "r")

  --Skip the first line (this is a title)
  c_stdio.fscanf(FIELD,  "%*[^\n]\n")

  finish = true
  while finish do
    var err = c_stdio.fgets(record, 200, FIELD)
    if err == [&int8](0) then
        finish = false
    end
    --Pull the first word
    key = io_utils.get_word(record, 1)
    --Convert to lowercase
    io_utils.to_lowercase(&key)
    if c_string.strncmp(key, "close", 5) == 0 then
        finish = false
        break
    elseif c_string.strncmp(key, "species", 7)  == 0 then
        for i=0, config[0].nspe do
            c_stdio.fgets(record1, 200, FIELD)
            regentlib.c.free(word)
            word = io_utils.get_word(record1, 1)
            if c_string.strlen(word) > 8 then
                regentlib.assert(false, "Name of species too long (max length is 8 characters)")
            end
            c_string.strcpy(config[0].namspe[i], word)
            config[0].masstmp[i] = io_utils.get_double(record1, 2)
            config[0].chgetmp[i] = io_utils.get_double(record1, 3)
            --TODO: NYI FROZEN Ignoring frozen for now
        end
    elseif c_string.strncmp(key, "bond", 4) == 0 then
        --TODO: NYI
        regentlib.assert(false, "bond keyword in FIELD file not yet supported")
    elseif c_string.strncmp(key, "cons", 4) == 0 then
        --TODO: NYI
        regentlib.assert(false, "cons keyword in FIELD file not yet supported")
    elseif c_string.strncmp(key, "angle", 5) == 0 then
        --TODO: NYI
        regentlib.assert(false, "angle keyword in FIELD file not yet supported")
    elseif c_string.strncmp(key, "dihed", 5) == 0 then
        --TODO: NYI
        regentlib.assert(false, "dihed keyword in FIELD file not yet supported")
    end 
  end


  c_stdio.fclose(FIELD)

--  @mxprm = imxprm
--  @cutoff = rcut
--  @nspe = inspe
--
--  @namspe = temp_namspe
--  @masstmp = temp_masstmp
--  @chgetmp = temp_chgetmp

  regentlib.c.free(word)
end

local fabs = regentlib.fabs(double)

task dl_meso_read_mod.read_field --( nspe : int, masstmp : &double, chgetmp : &double, npot : int, mxprm : int, namspe : &&int8,
                -- eunit : double, gamma : &double ) --, ...)
             ( config : region(ispace(int1d), config_type)) where writes(config), reads(config) do

  var FIELD = c_stdio.fopen("FIELD", "r")
  --Skip the first line (this is a title)
  c_stdio.fscanf(FIELD,  "%*[^\n]\n")
    --TODO NYI: engunit
    --if config[0].engunit == 1 then
    --    config[0].eunit = config[0].temp
    --else
        config[0].eunit = 1.0
    --end
  var record : int8[201]
  var record1 : int8[201]
  var key : &int8 = [&int8](regentlib.c.malloc(64))
  var word : &int8 = [&int8](regentlib.c.malloc(1))
  var word1 : &int8 = [&int8](regentlib.c.malloc(1))
  var spename : &int8 = [&int8](regentlib.c.malloc(9))
  var finmol : int

  var ispe: int
  var jspe : int
  var k : int
  var i : int
  var j : int
 
  var aa : double
  var bb : double
  var cc : double
  var dd : double
  var ee : double
  var ff : double
  var gg : double

  var interact : &&bool
  interact = [&&bool](regentlib.c.malloc([terralib.sizeof( &bool )] * 3))
  interact[0] = [&bool](regentlib.c.malloc([terralib.sizeof(bool)] * config[0].npot))
  interact[1] = [&bool](regentlib.c.malloc([terralib.sizeof(bool)] * config[0].npot))
  interact[2] = [&bool](regentlib.c.malloc([terralib.sizeof(bool)] * config[0].npot))
--Iniitalise interact to false
    for i=0, config[0].npot do
        interact[0][i] = false
        interact[1][i] = false
        interact[2][i] = false
    end
--  var vvv : &&double
--  vvv = [&&double](regentlib.c.malloc([terralib.sizeof(&double)] * mxprm))
--  for i=0, mxprm do
--    vvv[i] = [&double](regentlib.c.malloc([terralib.sizeof(double)] * npot))
--  end
--
--  var ktype : &int
--  ktype = [&int](regentlib.c.malloc([terralib.sizeof(int)] * npot))

  var finish : bool = true

  while finish do

    --Read line
    var err = c_stdio.fgets(record, 200, FIELD)
    if err == [&int8](0) then
        finish = false
    end
    regentlib.c.free(key)
    --Pull the first word
    key = io_utils.get_word(record, 1)
    --Convert to lowercase
    io_utils.to_lowercase(&key)
    if c_string.strncmp(key, "close", 5) == 0 then
        finish = false
        regentlib.c.free(key)
        break
    elseif c_string.strncmp(key, "species", 7) == 0 then
        for i = 0, config[0].nspe do
            c_stdio.fgets(record1, 200, FIELD)
            regentlib.c.free(word)
            word = io_utils.get_word(record1, 4)
            config[0].nspec[i] = c_stdlib.atoi(word)
            config[0].nusystcell = config[0].nusystcell + config[0].nspec[i]
        end 
        config[0].nsystcell = config[0].nsystcell + config[0].nusystcell
    elseif c_string.strncmp(key, "molecul", 7) == 0 then
      --TODO: NYI Molecules
        regentlib.assert(false, "molecul keyword in FIELD file not yet supported")
    elseif c_string.strncmp(key, "interact", 8) == 0 then
        regentlib.c.free(word)
        word = io_utils.get_word(record, 2)
        finmol = c_stdlib.atoi(word)
        for i = 0, finmol do
            err = c_stdio.fgets(record1, 200, FIELD)
            regentlib.c.free(word1)
            word1 = io_utils.get_word(record1, 1)
            c_string.strncpy( spename, word1, 8)
            spename[8] = int8(0)
            ispe = -1
            for j=0,config[0].nspe do
                if c_string.strcmp(spename, config[0].namspe[j]) == 0 then
                    ispe = j
                end
            end
            if ispe == -1 then
                regentlib.c.free(word1)
                word1 = io_utils.get_word(record1, 1)
                ispe = c_stdlib.atoi(word1)
            end
            if not (ispe >= 0 and ispe < config[0].nspe) then
                format.println("Failed to read FIELD error 47")
                c_stdlib.exit(47) 
            end

            regentlib.c.free(word1)
            word1 = io_utils.get_word(record1, 2)
            c_string.strncpy( spename, word1, 8)
            spename[8] = int8(0)
            jspe = -1
            for j=0, config[0].nspe do
                if c_string.strcmp(spename, config[0].namspe[j]) == 0 then
                    jspe = j
                end
            end
            if jspe == -1 then
                regentlib.c.free(word1)
                word1 =  io_utils.get_word(record1, 2)
                jspe = c_stdlib.atoi(word1)
            end
            if not (jspe >= 0 and jspe < config[0].nspe) then
                format.println("Failed to read FIELD error 47")
                c_stdlib.exit(47) 
            end

            --FIXME: This is a bit weird due to conversion to C (0-based) from Fortran - need to check. Looks good.
            if (ispe>jspe) then
                k = (( (ispe+1) * ispe) / 2 + jspe+1) - 1
            else
                k = (( (jspe+1) * jspe) / 2 + ispe+1) - 1
            end

            --TODO: NYI - handle frozen atoms, line 2373 in dl_meso


            aa = io_utils.get_double(record1, 4)
            bb = io_utils.get_double(record1, 5)
            cc = io_utils.get_double(record1, 6)
            dd = io_utils.get_double(record1, 7)
            ee = io_utils.get_double(record1, 8)
            ff = io_utils.get_double(record1, 9)
            gg = io_utils.get_double(record1, 10)

            regentlib.c.free(word1)
            word1 =  io_utils.get_word(record1, 3)
            io_utils.to_lowercase(&word1)
            if c_string.strncmp(word1, "lj", 2) == 0 then
                 --TODO: NYI
                regentlib.assert(false, "lj interaction in FIELD file not yet supported")
            elseif c_string.strncmp(word1, "wca", 3) == 0 then
                --TODO: NYI
                regentlib.assert(false, "wca interaction in FIELD file not yet supported")
            elseif c_string.strncmp(word1, "dpd", 3) == 0 then
                config[0].ktype[k] = 2
                config[0].vvv[0][k] = config[0].eunit * fabs(aa)
                config[0].vvv[1][k] = fabs(bb)
                config[0].vvv[2][k] = bb * bb
                config[0].gamma[k] = fabs(cc)
                interact[0][k] = true
                if fabs(bb) > 0.0 then
                    interact[1][k] = true
                end
                interact[2][k] = true
            elseif c_string.strncmp(word1, "mdpd", 4) == 0 then
                --TODO: NYI
                regentlib.assert(false, "mdpd interaction in FIELD file not yet supported")
            elseif c_string.strncmp(word1, "tab", 3) == 0 then
                --TODO: NYI
                regentlib.assert(false, "tab interaction in FIELD file not yet supported")
            end
            
        end -- end of interact
    elseif c_string.strncmp(key, "froz", 4) == 0 then
        --TODO: NYI Frozen atoms
        regentlib.assert(false, "froz keyword in FIELD file not yet supported")
    elseif c_string.strncmp(key, "surf", 4) == 0 then
        --TODO: NYI Surfaces
        regentlib.assert(false, "surf keyword in FIELD file not yet supported")
    elseif c_string.strncmp(key, "extern", 6) == 0 then
        --TODO: NYI external
        regentlib.assert(false, "extern keyword in FIELD file not yet supported")
    end
  end
  c_stdio.fclose(FIELD)
  --TODO: NYI Read tabulated interactions

  --Check for missing interactions and fix
  if (io_utils.ANY(interact, 0, config[0].npot, false)) then
    --TODO: NYI Check if lmb is enforced (line 2561)
    if false then
        --TODO: NYI
        format.println("Failed to read FIELD error 50")
        c_stdlib.exit(50) 
    else
        for i=0, config[0].nspe-1 do
            for j=i+1, config[0].nspe do
                ispe = (( (i+1) * (i+2)) / 2) - 1
                if( (not interact[0][ispe]) and ( not interact[1][ispe]) and ( not interact[2][ispe] )) then
                    format.println("Failed to read FIELD error 51")
                    c_stdlib.exit(51) 
                end
    
                jspe = (( (j+1) * (j+2)) / 2) - 1
                if( (not interact[0][jspe]) and ( not interact[1][jspe]) and ( not interact[2][jspe] )) then
                    format.println("Failed to read FIELD error 51")
                    c_stdlib.exit(51) 
                end
                k = (( (j+1) * j) /2 + (i+1)) - 1               
                --TODO: NYI Frozen (line 2582)
                if( (not interact[0][k]) or (not interact[1][k]) or (not interact[2][k])) then 
                    --Continue from line 2588
                    var typ = regentlib.fmax(config[0].ktype[ispe], config[0].ktype[jspe] )
                    config[0].ktype[k] = typ
                    if typ == -1 then
                        format.println("Failed to read FIELD error 68")
                        c_stdlib.exit(68) 
                    elseif typ == 0 then
                        --TODO: NYI
                        format.println("lj not supported")
                    elseif typ == 1 then
                        --TODO: NYI
                        format.println("typ 1 not supported")
                    elseif typ == 2 then
                        if interact[0][k] then 
                            aa = config[0].vvv[0][k]
                        else
                            aa = regentlib.c.sqrt( config[0].vvv[0][ispe] * config[0].vvv[0][jspe] )
                            config[0].vvv[0][k] = aa
                        end
                        if interact[1][k] then
                            bb = config[0].vvv[1][k]
                        else
                            bb = 0.5 * ( config[0].vvv[1][ispe] + config[0].vvv[1][jspe] )
                            config[0].vvv[1][k] = bb
                            config[0].vvv[2][k] = bb * bb
                        end
                        if( not interact[2][k] ) then
                            config[0].gamma[k] = regentlib.c.sqrt( config[0].gamma[ispe] * config[0].gamma[jspe] )
                        end
                        interact[0][k] = true
                        interact[1][k] = true
                        interact[2][k] = true
                    end
                end
    
                               
 
            end
        end
    end
  end
  regentlib.c.free(interact[0])
  regentlib.c.free(interact[1])
  regentlib.c.free(interact[2])
  regentlib.c.free(interact)
  --No need to broadcast results - we will use Regions to do this.

--TODO: Support nfold
--!     determine numbers of particles, bonds etc. from unit cell values
--
      config[0].nsyst = config[0].nsystcell -- * nfold
      config[0].nusyst = config[0].nusystcell -- * nfold
--      nfsyst = nfsystcell * nfold
--      numbond = numbondcell * nfold
--      numcon = numconcell * nfold
--      numang = numangcell * nfold
--      numdhd = numdhdcell * nfold
--      numexc = numexccell * nfold
        config[0].nspec = config[0].nspec --* nfold
--      nspecmol = nspecmol * nfold
--      nummol = nummolcell * nfold

--TODO: NYI Rescale bond energy parameters

--TODO: NYI - rescale eunit
--!     rescale electrostatic permittivity parameter
--
--      gammaelec = gammaelec * eunit
--      bjerelec = bjerelec * eunit

  regentlib.c.free(word)
end

task dl_meso_read_mod.scan_control( config : region(ispace(int1d), config_type)) where writes(config.l_exist, config.l_safe, config.l_scr,
                    config.l_temp, config.l_time, config.l_conf, config.l_init, config.l_rest, config.temp, config.tstep) 
do

    config[0].l_scr = false
    config[0].l_temp = false
    config[0].l_time = false
    config[0].l_safe = false
    config[0].l_conf = true
    config[0].l_init = false
    config[0].l_rest = false
    var record : int8[201]
    var key : &int8 = [&int8](regentlib.c.malloc(64))
    var key1 : &int8 = [&int8](regentlib.c.malloc(64))

    if( c_unistd.access("CONTROL", c_unistd.R_OK) ~= 0 ) then
        config[0].l_exist = false
    else
        config[0].l_exist = true
    end
    var CONTROL = c_stdio.fopen("CONTROL", "r")
    if CONTROL ~= [&c_stdio.FILE](int64(0)) then
        config[0].l_safe = true
        c_stdio.fscanf(CONTROL,  "%*[^\n]\n")
        var finish = true
        while finish do
            var err = c_stdio.fgets(record, 200, CONTROL)
            if err == [&int8](0) then
                finish = false
                break
            end
            regentlib.c.free(key)
            key = io_utils.get_word(record, 1)
            io_utils.to_lowercase(&key)

            if c_string.strncmp(key, "finish", 6) == 0 then 
            elseif c_string.strncmp(key, "l_scr", 5) == 0 then
                config[0].l_scr = true
            elseif c_string.strncmp(key, "l_init", 6) == 0 then
                config[0].l_init = true
            elseif c_string.strncmp(key, "no", 2) == 0 then
                regentlib.c.free(key1)
                io_utils.to_lowercase(&key1)
                key1[4] = int8(0)
                if c_string.strncmp(key1, "conf", 4) == 0 then
                    config[0].l_conf = false
                end
            elseif c_string.strncmp(key, "restart", 7) == 0 then
                config[0].l_rest = true
            elseif c_string.strncmp(key, "temp", 4) == 0 then
                var temp = io_utils.get_double(record, 2)
                config[0].temp = temp
                if temp > 1e-16 then
                    config[0].l_temp = true
                end
            elseif c_string.strncmp(key, "timestep", 8) == 0 then
                var tstep = io_utils.get_double(record, 2)
                config[0].tstep = tstep
                if tstep > 1e-16 then
                    config[0].l_time = true
                end
            end
        end
    end
    c_stdio.fclose(CONTROL)
    regentlib.c.free(key)
    regentlib.c.free(key1)
end

task dl_meso_read_mod.read_control(config : region(ispace(int1d), config_type), l_readvol : bool, l_config : bool, l_rest : bool) 
                    where reads(config), writes(config.tclose, config.cutoff,
                                                config.itype, config.btype, config.nseql, config.ldyn, config.nsbpo, config.nsbts,
                                                config.ltemp, config.nstk, config.iscorr, config.lcorr, config.nrun, config.temp,
                                                config.straj, config.ntraj, config.keytrj, config.ltraj, config.tstep, config.rtstep,
                                                config.tstepsq, config.volm, config.space, config.text)
do
    
    var record : int8[201]
    var key1 : &int8 = [&int8](regentlib.c.malloc(64))
    var key2 : &int8 = [&int8](regentlib.c.malloc(64))
    var word1 : &int8 = [&int8](regentlib.c.malloc(64))
    var word2 : &int8 = [&int8](regentlib.c.malloc(64))

    var safe : bool = true
    var lelec : bool = false
    var lspme : bool = false
    var lewprc : bool = false
    var lrpszero : bool = false
    var slbetlen : int = 0

    var CONTROL = c_stdio.fopen("CONTROL", "r")
    if CONTROL == [&c_stdio.FILE](int64(0)) then
        format.println("Failed to read control Error 21")
        c_stdlib.exit(21) 
    end
    --Skip the title line
--    c_stdio.fscanf(CONTROL,  "%*[^\n]\n")
    --Read the title line into the text field.
    c_stdio.fgets(record, 80, CONTROL)
    --Can't call strcpy on fields of the region (no aliasing allowed)
    for i = 0, 81 do
        config[0].text[i] = record[i]
    end
    var finish : bool = true
    while finish do
        var err = c_stdio.fgets(record, 200, CONTROL)
        if err == [&int8](0) then
            finish = false
            break
        end
        regentlib.c.free(key1)
        regentlib.c.free(key2)
        key1 = io_utils.get_word(record, 1)
        key2 = io_utils.get_word(record, 2)
        io_utils.to_lowercase(&key1)
        io_utils.to_lowercase(&key2)
        if c_string.strncmp(key1, "finish", 6) == 0 then
            finish = false
            break
        elseif c_string.strncmp(key1, "bjer", 4) == 0 then
            --TODO: NYI bjer...
            regentlib.assert(false, "bjer keyword not yet supported.")
        elseif c_string.strncmp(key1, "bound", 5) == 0 and c_string.strncmp(key2, "halo", 4) == 0 then
            --TODO NYI: Boundary halo
            regentlib.assert(false, "bound keyword not yet supported")
        elseif c_string.strncmp(key1, "conf", 4) == 0 then
            --TODO NYI: conf...
            regentlib.assert(false, "conf keyword not yet supported")
        elseif c_string.strncmp(key1, "close", 5) == 0 and c_string.strncmp(key2, "time", 4) == 0 then
            config[0].tclose = io_utils.get_double(record, 3)
        elseif c_string.strncmp(key1, "cut", 3) == 0 or c_string.strncmp(key2, "rcut", 4 ) == 0 then
            config[0].cutoff = io_utils.get_double(record, 2)
        elseif c_string.strncmp(key1, "densvar", 7) == 0 then
            --TODO NYI: densvar
            regentlib.assert(false, "densvar keyword not yet supported")
        elseif c_string.strncmp(key1, "elec", 4) == 0 and c_string.strncmp(key2, "cut", 3) == 0 then
            --TODO NYI: elec... cut...
            regentlib.assert(false, "elec keyword not yet supported")
        elseif c_string.strncmp(key1, "ensemble", 8) == 0 then
            regentlib.c.free(word1)
            word1 = io_utils.get_word(record, 3)
            io_utils.to_lowercase(&word1)
            var compword : int8[5]
            compword[0] = word1[0]
            compword[1] = word1[1]
            compword[2] = word1[2]
            compword[3] = word1[3]
            compword[4] = int8(0)
            if c_string.strcmp(compword, "mdvv") == 0 then
                config[0].itype = 0
                regentlib.c.free(word2)
                word2 = io_utils.get_word(record, 4)
            elseif c_string.strcmp(compword, "dpdv") == 0 then
                --TODO NYI: dpdv
                regentlib.assert(false, "dpdv not yet supported")
            elseif c_string.strcmp(compword, "dpds") == 0 then
                --TODO NYI: dpds
                regentlib.assert(false, "dpds not yet supported")
            elseif c_string.strcmp(compword, "lowe") == 0 then
                --TODO NYI: lowe
                regentlib.assert(false, "lowe not yet supported")
            elseif c_string.strcmp(compword, "pete") == 0 then
                --TODO NYI: pete
                regentlib.assert(false, "pete not yet supported")
            elseif c_string.strcmp(compword, "stoy") == 0 then
                regentlib.assert(false, "stoy not yet supported")
                --TODO NYI: stoy
            end
            --Continue from line 1180
            compword[0] = key2[0]
            compword[1] = key2[1]
            compword[2] = key2[2]
            compword[3] = int8(0)
            if c_string.strcmp(compword, "nvt") == 0 then
                config[0].btype = 0
            elseif c_string.strcmp(compword, "npt") == 0 then 
                --TODO NYI: npt
                regentlib.assert(false, "npt not yet supported")
            elseif c_string.strcmp(compword, "nst") == 0 then
                --TODO NYI: nst
                regentlib.assert(false, "nst not yet supported")
            end 
        elseif c_string.strncmp(key1, "equil", 5) == 0 then
            if c_string.strncmp(key2, "steps", 5) == 0 then
                regentlib.c.free(word1)
                word1 = io_utils.get_word(record, 3)
                config[0].nseql = c_stdlib.atoi(word1)
            else
                regentlib.c.free(word1)
                word1 = io_utils.get_word(record, 2)
                config[0].nseql = c_stdlib.atoi(word1)
            end
        elseif c_string.strncmp(key1, "ewald", 5) == 0 then
            --TODO NYI: ewald
                regentlib.assert(false, "ewald not yet supported")
        elseif c_string.strncmp(key1, "froz", 4) == 0 then
                regentlib.assert(false, "froz not yet supported")
            --TODO NYI: froz
        elseif c_string.strncmp(key1, "global", 6) == 0 and c_string.strncmp(key2, "bonds", 5) == 0 then
            --TODO NYI: global bonds
                regentlib.assert(false, "global bonds not yet supported")
        elseif c_string.strncmp(key1, "io", 2) == 0 and c_string.strncmp(key2, "writ", 4) == 0 then
            --TODO NYI: io write
                regentlib.assert(false, "io write not yet supported")
        elseif c_string.strncmp(key1, "job", 3) == 0 and c_string.strncmp(key2, "time", 4) == 0 then
            --TODO: We ignore job time for now
            --config[0].timjob = io_utils.get_double(record, 3)
            format.println("WARNING: Ignoring job time from CONTROL file at this time")
        elseif c_string.strncmp(key1, "many", 4) == 0 and c_string.strncmp(key2, "cut", 3) == 0 then
            --TODO NYI: many cut
                regentlib.assert(false, "many cut not yet supported")
        elseif c_string.strncmp(key1, "mxshak", 6) == 0 then
            --TODO NYI: mxshak
                regentlib.assert(false, "mxshak not yet supported")
        elseif c_string.strncmp(key1, "ndump", 5) == 0 then
            --TODO NYI: ndump
                regentlib.assert(false, "ndump not yet supported")
        elseif c_string.strncmp(key1, "nfold", 5) == 0 and (l_config or l_rest) then
            --TODO NYI: nfold
                regentlib.assert(false, "nfold not yet supported")
        elseif c_string.strncmp(key1, "no", 2) == 0 then
            --TODO NYI: no
                regentlib.assert(false, "no keyword options not yet supported")
        elseif c_string.strncmp(key1, "openmp", 6) == 0 then
            --Not supporting openmp flag - parallelism already handled
            format.println("WARNING: Ignoring openmp request, inbuilt parallelism used instead")
        elseif c_string.strncmp(key1, "perm", 4) == 0 then
            --TODO NYI: perm
                regentlib.assert(false, "perm not yet supported")
        elseif c_string.strncmp(key1, "pres", 4) == 0 then
            --TODO NYI: pers
                regentlib.assert(false, "pres not yet supported")
        elseif c_string.strncmp(key1, "print", 5) == 0 then
            if c_string.strncmp(key2, "part", 4) == 0 then
                config[0].ldyn = true
            else
                if c_string.strncmp(key2, "every", 5) == 0 then
                    regentlib.c.free(word1)
                    word1 = io_utils.get_word(record, 3)
                    config[0].nsbpo = c_stdlib.atoi(word1)
                else
                    regentlib.c.free(word1)
                    word1 = io_utils.get_word(record, 2)
                    config[0].nsbpo = c_stdlib.atoi(word1)
                end
            end
        elseif c_string.strncmp(key1, "restart", 7) == 0 then
            --TODO NYI: restart
                regentlib.assert(false, "restart not yet supported")
        elseif c_string.strncmp(key1, "scale", 5) == 0 then
            regentlib.c.free(word1)
            word1 = io_utils.get_word(record, 2)
            config[0].nsbts = c_stdlib.atoi(word1)
            if config[0].nsbts == 0 then
                regentlib.c.free(word1)
                word1 = io_utils.get_word(record, 3)
                config[0].nsbts = c_stdlib.atoi(word1)
            end
            if config[0].nsbts == 0 then
                regentlib.c.free(word1)
                word1 = io_utils.get_word(record, 4)
                config[0].nsbts = c_stdlib.atoi(word1)
            end
            if config[0].nsbts > 0 then
                config[0].ltemp = true
            end
        elseif c_string.strncmp(key1, "seed", 4) == 0 then
                regentlib.assert(false, "seed not yet supported")
            --TODO NYI: seed
        elseif c_string.strncmp(key1, "shake", 5) == 0 then
                regentlib.assert(false, "shake not yet supported")
            --TODO NYI: shake
        elseif c_string.strncmp(key1, "smear", 5) == 0 then
            --TODO NYI: smear
                regentlib.assert(false, "smear not yet supported")
        elseif c_string.strncmp(key1, "spme", 4) == 0 then
            --TODO NYI: spme
                regentlib.assert(false, "spme not yet supported")
        elseif c_string.strncmp(key1, "stack", 5) == 0 then
            if c_string.strncmp(key2, "size", 4) == 0 then
                regentlib.c.free(word1)
                word1 = io_utils.get_word(record, 3)
                config[0].nstk = c_stdlib.atoi(word1)
            else
                regentlib.c.free(word1)
                word1 = io_utils.get_word(record, 2)
                config[0].nstk = c_stdlib.atoi(word1)
            end
            regentlib.assert(config[0].nstk < CONST_maxstk, "The Stack size is larger than the maximum - please increase CONST_maxstk to continue")
        elseif c_string.strncmp(key1, "stats", 5) == 0 then
            if c_string.strncmp(key2, "every", 5) == 0 then
                regentlib.c.free(word1)
                word1 = io_utils.get_word(record, 3)
                config[0].iscorr = c_stdlib.atoi(word1)
            else
                regentlib.c.free(word1)
                word1 = io_utils.get_word(record, 2)
                config[0].iscorr = c_stdlib.atoi(word1)
            end
            if config[0].iscorr > 0 then
                config[0].lcorr = true
            end
        elseif c_string.strncmp(key1, "steps", 5) == 0 then
            config[0].nrun = c_stdlib.atoi(key2)
        elseif c_string.strncmp(key1, "stres", 5) == 0 then
            --TODO NYI: stres
                regentlib.assert(false, "stres not yet supported")
        elseif c_string.strncmp(key1, "surf", 4) == 0 then
            --TODO NYI: surf
                regentlib.assert(false, "surf not yet supported")
        elseif c_string.strncmp(key1, "temp", 4) == 0 then
            config[0].temp = io_utils.get_double(record, 2)
        elseif c_string.strncmp(key1, "therm", 5) == 0 and c_string.strncmp(key2, "cut", 3) == 0 then
            --TODO NYI: therm cut
                regentlib.assert(false, "therm cut not yet supported")
        elseif c_string.strncmp(key1, "traj", 4) == 0 then
            regentlib.c.free(word1)
            word1 = io_utils.get_word(record, 2)
            config[0].straj = c_stdlib.atoi(word1)        
            regentlib.c.free(word1)
            word1 = io_utils.get_word(record, 3)
            config[0].ntraj = c_stdlib.atoi(word1)        
            regentlib.c.free(word1)
            word1 = io_utils.get_word(record, 4)
            config[0].keytrj = c_stdlib.atoi(word1)
            if config[0].straj == 0 then
                config[0].straj = config[0].nseql
            end
            if config[0].keytrj < 0 then 
                config[0].keytrj = 0
            end
            if config[0].keytrj > 2 then
                config[0].keytrj = 2
            end
            if config[0].ntraj > 0 then
                config[0].ltraj = true
            end
        elseif c_string.strncmp(key1, "timestep", 8) == 0 then
            config[0].tstep = io_utils.get_double(record, 2)
        elseif c_string.strncmp(key1, "vacu", 4) == 0 then
            --TODO NYI: vacu
            regentlib.assert(false, "vacu keyword in CONTROL file not yet supported")
        elseif c_string.strncmp(key1, "vol", 3) == 0 and l_readvol then
            config[0].space.dim_x = io_utils.get_double(record, 2)    
            config[0].space.dim_y = io_utils.get_double(record, 3)    
            config[0].space.dim_z = io_utils.get_double(record, 4)
            if config[0].space.dim_x * config[0].space.dim_y * config[0].space.dim_z < 1e-10 then
                var cubeside = c_math.cbrtf(config[0].space.dim_x)
                config[0].space.dim_x = cubeside            
                config[0].space.dim_y = cubeside            
                config[0].space.dim_z = cubeside            
            end   
            if not l_config then
                --TODO NYI: lnfold
            else
                --TODO NYI: lnfold
            end 
        end
    end

    c_stdio.fclose(CONTROL)

--!     check that volume has been read either in CONTROL, CONFIG or export file:
--!     if using CONFIG file but not restarting, set total size of system according
--!     to nfold duplication
    if l_config and (not l_rest) then
        if config[0].space.dim_x * config[0].space.dim_y * config[0].space.dim_z < 1.0e-10 then
            format.println("Failed to read config file: ERROR 9")
            c_stdlib.exit(9)
        end
        --TODO NYI: lnfold
    else
        if config[0].space.dim_x * config[0].space.dim_y * config[0].space.dim_z < 1.0e-10 then
            format.println("Failed to read config file: ERROR 9")
            c_stdlib.exit(9)
        end
    end

    --TODO NYI: electrostatic stuff - line 1625 to 1648
    --TODO NYI: beta value stuff - line 1650 to 1680
    --TODO NYI: ewald convergence stuff - line 1682 to 1723
    --TODO NYI: electrostatic type stuff - line 1725 to 1733
    --TODO NYI: Check pressure for ensembles with barostats - line 1735 to 1739

--Don't need to broadcast - region handles this

    --TODO NYI: Nfold set - line 1915
    config[0].rtstep = 1.0 / config[0].tstep
    config[0].tstepsq = config[0].tstep * config[0].tstep
    config[0].volm = config[0].space.dim_x * config[0].space.dim_y * config[0].space.dim_z
    --TODO NYI: frozen things - line 1919
    config[0].nsbpo = regentlib.fmax(config[0].nsbpo, 1)


    regentlib.c.free(key1)
    regentlib.c.free(key2)
    regentlib.c.free(word2)
    regentlib.c.free(word1)
end


return dl_meso_read_mod

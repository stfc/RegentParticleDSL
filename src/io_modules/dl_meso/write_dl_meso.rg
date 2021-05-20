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

task write_output_summary( time : double, lbegin : bool, l_scr : bool, config : region(ispace(int1d), config_type) ) where reads(config.stpte,
                            config.stppe, config.stpvir, config.stpprs, config.stpttp, config.rav, config.nstep, config.stptke)
do

    --We don't keep file open here
    var OUTPUT = c_stdio.fopen('OUTPUT', 'a')
    var i : int
    --Assumed for now, config[0].outsel is 0
    --if config[0].outsel == 0 then
        --Static system without bonds, angles, dihedrals, barostat, electrostatics, surfaces
        if lbegin then
            c_stdio.fprintf(OUTPUT, "\n ")
            for i=0, 95 do
                c_stdio.fprintf(OUTPUT, "-")
            end
            c_stdio.fprintf(OUTPUT, "\n       ")
            c_stdio.fprintf(OUTPUT, "step      en-total      pe-total     vir-total      ke-total      pressure   temperature\n ")
            for i=0, 95 do
                c_stdio.fprintf(OUTPUT, "-")
            end
            c_stdio.fprintf(OUTPUT, "\n")
        end
        var buffer : int8[25]
        c_stdio.snprintf(&(buffer[0]), 25, "%i", config[0].nstep)
        c_stdio.fprintf(OUTPUT, " ")
        c_stdio.fprintf(OUTPUT, "%10s", buffer)
        c_stdio.fprintf(OUTPUT, "  ")
        c_stdio.fprintf(OUTPUT, "%12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", config[0].stpte, config[0].stppe,
                              config[0].stpvir, config[0].stptke, config[0].stpprs, config[0].stpttp )
        c_stdio.fprintf(OUTPUT, " %10.3f  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", time,
                              config[0].rav[0], config[0].rav[1], config[0].rav[7], config[0].rav[8],
                              config[0].rav[9], config[0].rav[12] )
        c_stdio.fprintf(OUTPUT, " ")
        for i=0, 95 do
            c_stdio.fprintf(OUTPUT, "-")
        end
        c_stdio.fprintf(OUTPUT, "\n")
    --elseif ....
    --end
    if not l_scr then   
        --TODO NYI: Don't quite understand this code (nor how to do it in C)
        --Think since we're closing the file its not needed.
        --TODO NYI: Don't yet support l_scr
    end
    c_stdio.fclose(OUTPUT)
end


task write_output_equil( config : region(ispace(int1d), config_type) ) where reads(config.nstep) 
do

    var OUTPUT = c_stdio.fopen('OUTPUT', 'a')
    var i : int = 0
    --if config[0].outsel == 0 then
        var buffer : int8[25]
        c_stdio.snprintf(&(buffer[0]), 25, "%i", config[0].nstep)
        c_stdio.fprintf(OUTPUT, "\n equilibration period ended at step %10s\n\n ", buffer) 
        for i=0,95 do
            c_stdio.fprintf(OUTPUT, "-")
        end
        c_stdio.fprintf(OUTPUT, "\n")
    --elseif ....
    --end
    c_stdio.fclose(OUTPUT)
end

task write_correl( config : region(ispace(int1d), config_type) , time : double) where reads(config.CORREL_newjob, config.stress,
                    config.stpte, config.stppe, config.stpprs, config.stpttp),  writes(config.CORREL_newjob)
do
    if(config[0].CORREL_newjob) then
        var exists = false
  --Check FIELD file exists
        if( c_unistd.access("CORREL", c_unistd.W_OK) ~= 0 ) then
            exists = true
        end
        var CORREL : &c_stdio.FILE
        if exists then
            CORREL = c_stdio.fopen('CORREL', 'a')
        else
            CORREL = c_stdio.fopen('CORREL', 'w')
            --Write header for CORREL file
            --if config[0].outsel == 0 then
                c_stdio.fprintf(CORREL, "#         time      en-total      pe-total      pressure          s_xx          s_xy          s_xz          s_yx          s_yy          s_yz          s_zx          s_z    y          s_zz   temperature\n")
            --elseif ....
            --end
        end
        config[0].CORREL_newjob = false
        c_stdio.fclose(CORREL)
    end
    var totstress : double[9]
    var i : int
    for i=0, 9 do
        var lowindex = 4*(i+1) - 3 -1
        var hiindex = 4*(i+1)
        var j : int
        totstress[i] = 0.0
        for j = lowindex, hiindex do
            totstress[i] = totstress[i] + config[0].stress[j]
        end
    end
    --if config[0].outsel == 0 then
       var CORREL = c_stdio.fopen('CORREL', 'a')  
        c_stdio.fprintf(CORREL,  "   %12.6e", time)
        c_stdio.fprintf(CORREL, "  %12.6e", config[0].stpte)
        c_stdio.fprintf(CORREL, "  %12.6e", config[0].stppe)
        c_stdio.fprintf(CORREL, "  %12.6e", config[0].stpprs)
        c_stdio.fprintf(CORREL, "  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e", totstress[0],
                        totstress[1], totstress[2], totstress[3], totstress[4], totstress[5], totstress[6], totstress[7], totstress[8])
        c_stdio.fprintf(CORREL, "  %12.6e\n", config[0].stpttp)
    --elseif ....
    --end
    c_stdio.fclose(CORREL)
end

task write_stress(time : double)
  regentlib.assert(false, "Can't yet write stress files")
end

--TODO: Assuming not needed but not sure
--task init_output_groups()
--end

task write_history_header(config : region(ispace(int1d), config_type), parts : region(ispace(int1d), part)) where 
                                                                                reads(config.nstep, config.kres, config.text,
                                                                                    config.nspe, config.nusyst, config.nsyst,
                                                                                    config.ktype, config.namspe, config.filesize,
                                                                                    config.numframe, config.masstmp, config.chgetmp,
                                                                                    config.vvv, config.markerpos, config.headersize),
                                                                                reads(parts.lab, parts.ltp, parts.ltm),
                                                                                writes(config.headersize, config.numframe, config.markerpos, config.filesize)
do

    var lhist : bool = c_unistd.access("HISTORY", c_unistd.W_OK) ~= 0
    var nstep1 : int = 0
    var newfile : bool = (config[0].nstep==0 or config[0].kres>1  or (not lhist))
--    newfile = false
    if not newfile then
        var HISTORY = c_stdio.fopen('HISTORY', "r+b")
        var buffer : int32[32]
        c_stdio.fread(&(buffer[0]),[terralib.sizeof(int32)],2, HISTORY)
        --TODO NYI: Other restart related-things

        c_stdio.fclose(HISTORY)
    else
        --Starting new HISTORY file
        config[0].filesize = 0
        config[0].numframe = 0
        nstep1 = 0
        var HISTORY = c_stdio.fopen('HISTORY', 'w+b')
        --For now we're forcing endianness is 1 and double precision
        var buffer : int32[32]
        buffer[0] = 1
        buffer[1] = 8
        c_stdio.fwrite(&(buffer[0]), [terralib.sizeof(int32)], 2, HISTORY);
        @([&int64](&(buffer[0]))) = config[0].filesize
--        ([&int64](buffer[0]))[0] = config[0].filesize
        buffer[2] = config[0].numframe
        buffer[3] = nstep1
        c_stdio.fwrite(&(buffer[0]), [terralib.sizeof(int32)], 4, HISTORY)
        --Write the name of the simulation
        var name : int8[81]
        c_string.strcpy(&(name[0]), config[0].text)
        c_stdio.fwrite( &(name[0]), [terralib.sizeof(int8)], 80, HISTORY)
        --write numbers of species, molecule types, unbonded beads, total
        --number of beads, number of bonds and constraints, data-writing
        --level and boundary types
        buffer[0] = config[0].nspe
        buffer[1] = 0
        buffer[2] = config[0].nusyst
        buffer[3] = config[0].nsyst
        buffer[4] = 0
        buffer[5] = 0
        buffer[6] = 0
        buffer[7] = 0
        buffer[8] = 0
        c_stdio.fwrite(&(buffer[0]), [terralib.sizeof(int32)], 9, HISTORY)
        
        --write names, masses, interaction lengths, charges and frozen property for all species
        --For Fortran compatibility we always ignore the 9th element (which is \0 if the string is 8 long)
        for i=0, config[0].nspe do
            var k = ( (i+1) *(i+1)) / 2 - 1
            if config[0].ktype[k] == -1 then
                --Write name
                c_string.strcpy(&(name[0]), config[0].namspe[i])
                c_stdio.fwrite(&(name[0]), [terralib.sizeof(int8)],8,HISTORY)
                --Write mass
                var mass = config[0].masstmp[i]
                c_stdio.fwrite(&mass, [terralib.sizeof(double)], 1, HISTORY)
                --Write interaction length?
                var vvv = config[0].vvv[0][k]
                c_stdio.fwrite(&vvv, [terralib.sizeof(double)], 1, HISTORY)
                --Write charge
                var chge = config[0].chgetmp[i]
                c_stdio.fwrite(&chge, [terralib.sizeof(double)], 1, HISTORY)
                --Never frozen yet TODO: NYI
                buffer[0] = 0
                c_stdio.fwrite(&(buffer[0]), [terralib.sizeof(int32)], 1, HISTORY)
            elseif config[0].ktype[k] < 3 then
                --Write name
                c_string.strcpy(&(name[0]), config[0].namspe[i])
                c_stdio.fwrite(&(name[0]), [terralib.sizeof(int8)], 8 ,HISTORY)
                --Write mass
                var mass = config[0].masstmp[i]
                c_stdio.fwrite(&(mass), [terralib.sizeof(double)], 1, HISTORY)
                --Write interaction length?
                var vvv = config[0].vvv[1][k]
                c_stdio.fwrite(&vvv, [terralib.sizeof(double)], 1, HISTORY)
                --Write charge
                var chge = config[0].chgetmp[i]
                c_stdio.fwrite(&chge, [terralib.sizeof(double)], 1, HISTORY)
                --Never frozen yet TODO: NYI
                buffer[0] = 0
                c_stdio.fwrite(&(buffer[0]), [terralib.sizeof(int32)], 1, HISTORY)
            elseif config[0].ktype[k] == 3 then
                --Write name
                c_string.strcpy(&(name[0]), config[0].namspe[i])
                c_stdio.fwrite(&(name[0]), [terralib.sizeof(int8)],8,HISTORY)
                --Write mass
                var mass = config[0].masstmp[i]
                c_stdio.fwrite(&mass, [terralib.sizeof(double)], 1, HISTORY)
                --Write interaction length?
                var vvv = config[0].vvv[5][k]
                c_stdio.fwrite(&vvv, [terralib.sizeof(double)], 1, HISTORY)
                --Write charge
                var chge = config[0].chgetmp[i]
                c_stdio.fwrite(&chge, [terralib.sizeof(double)], 1, HISTORY)
                --Never frozen yet TODO: NYI
                buffer[0] = 0
                c_stdio.fwrite(&(buffer[0]), [terralib.sizeof(int32)], 1, HISTORY)
            end
        end

        --TODO: NYI - don't write molcule names
        --Close the file
        c_stdio.fclose(HISTORY)

        --TODO NYI: Not yet supporting markerpos
        --int32 is 4 bytes
        var Ilen_li : int64 = 4
        var LIlen_li : int64 = 8
        var Dlen_li : int64 = 8
        config[0].markerpos = 2 * Ilen_li + 1
        var mypos : int64 = config[0].markerpos + 80 + LIlen_li + 11 + Ilen_li + config[0].nspe * ( 8 + Ilen_li + 3 * Dlen_li)
        
        var beadcount = config[0].nsyst
        --Collect bead information
        var collect_buffer : &int32 = [&int32](regentlib.c.malloc([terralib.sizeof(int32)] * beadcount*4))
        if collect_buffer == [&int32](int32(0)) then
            regentlib.assert(false, "Failed to write history file, error 1294")
        end
        for i = 0, beadcount do
            collect_buffer[4*i] = parts[i].lab
            collect_buffer[4*i +1] = parts[i].ltp
            collect_buffer[4*i +2] = parts[i].ltm
            collect_buffer[4*i +3] = 0
        end
        --Write this information to the history file
        HISTORY = c_stdio.fopen('HISTORY', 'r+b')
        --Move to the position expected in the history file
        c_stdio.fseek(HISTORY, mypos, c_stdio.SEEK_SET)
        --Write the collect_buffer array to the file
        c_stdio.fwrite(collect_buffer, [terralib.sizeof(int32)], 4*beadcount, HISTORY)
        --Close the file
        c_stdio.fclose(HISTORY)
        mypos = mypos + 4 * config[0].nsyst * Ilen_li 
        regentlib.c.free(collect_buffer)
        --TODO NYI: Write bonds and constraints

    end
    config[0].headersize = io_utils.get_file_size("HISTORY")
    config[0].filesize = config[0].headersize
    var HISTORY = c_stdio.fopen('HISTORY', 'r+b')
    c_stdio.fseek(HISTORY, config[0].markerpos, c_stdio.SEEK_SET)
    var hsize = config[0].headersize
    c_stdio.fwrite(&hsize, [terralib.sizeof(int64)], 1, HISTORY)
    var ti32 = config[0].numframe
    c_stdio.fwrite(&ti32, [terralib.sizeof(int32)], 1, HISTORY)
    ti32 =config[0].nstep
    c_stdio.fwrite(&ti32, [terralib.sizeof(int32)], 1, HISTORY)
    c_stdio.fclose(HISTORY)
end




task write_history(config : region(ispace(int1d), config_type), parts : region(ispace(int1d), part), time : double) where 
                                                                                reads(config),
                                                                                        reads(parts.core_part_space,
                                                                                        parts.lab, parts.fxx, parts.fyy, parts.fzz, parts.neighbour_part_space),
                                                                                writes(config.numframe, config.filesize)
do

    --int32 is 4 bytes
    var Ilen_li : int64 = 4
    var LIlen_li : int64 = 8
    var Dlen_li : int64 = 8
    var mypos = config[0].filesize
    config[0].filesize = config[0].filesize + int64(7)*Dlen_li + Ilen_li + int64(config[0].nsyst) * (Ilen_li + int64(3) * int64(config[0].keytrj+1)*Dlen_li)
    config[0].numframe = config[0].numframe + 1

    --Open HISTORY file, update header data and find beginning of new frame, and close file
    var HISTORY = c_stdio.fopen("HISTORY", "r+b") 
    c_stdio.fseek(HISTORY, config[0].markerpos, c_stdio.SEEK_SET)
    --Can't access pointers to fields...
    var ti64 = config[0].filesize
    c_stdio.fwrite(&ti64, [terralib.sizeof(int64)], 1, HISTORY)
    var ti32 = config[0].numframe
    c_stdio.fwrite(&ti32, [terralib.sizeof(int32)], 1, HISTORY)
    ti32 = config[0].nstep
    c_stdio.fwrite(&ti32, [terralib.sizeof(int32)], 1, HISTORY)

    c_stdio.fseek(HISTORY, mypos, c_stdio.SEEK_SET)
    c_stdio.fwrite(&time, [terralib.sizeof(double)], 1, HISTORY)
    var nsyst = config[0].nsyst
    c_stdio.fwrite(&nsyst, [terralib.sizeof(int32)], 1, HISTORY)
    var buf : double[6]
    buf[0] = config[0].space.dim_x
    buf[1] = config[0].space.dim_y
    buf[2] = config[0].space.dim_z
    buf[3] = config[0].shrdx
    buf[4] = config[0].shrdy
    buf[5] = config[0].shrdz
    c_stdio.fwrite(&(buf[0]), [terralib.sizeof(double)], 6, HISTORY)
   
    mypos = mypos + int64(7) * Dlen_li + Ilen_li 
    --Continue from 2478
    var double_buffer : &double = [&double](regentlib.c.malloc([terralib.sizeof(double)] * 3*(config[0].keytrj + 1) * config[0].nsyst))
    var labs : &int32 = [&int32](regentlib.c.malloc([terralib.sizeof(int32)] * config[0].nsyst))
    regentlib.assert(not isnull(double_buffer), "Failed to allocate double_buffer")
    regentlib.assert(not isnull(labs), "Failed to allocate labs")
    if config[0].keytrj == 0 then
        var ii = 0
        for i in parts.ispace do
            if neighbour_init.check_valid(parts[i].neighbour_part_space) then
                double_buffer[3*ii] = parts[i].core_part_space.pos_x
                double_buffer[3*ii+1] = parts[i].core_part_space.pos_y
                double_buffer[3*ii+2] = parts[i].core_part_space.pos_z
                labs[ii] = parts[i].lab
                ii = ii + 1
            end
        end
    elseif config[0].keytrj == 1 then
        var ii = 0
        for i in parts.ispace do
            if neighbour_init.check_valid(parts[i].neighbour_part_space) then
                double_buffer[6*ii] = parts[i].core_part_space.pos_x
                double_buffer[6*ii+1] = parts[i].core_part_space.pos_y
                double_buffer[6*ii+2] = parts[i].core_part_space.pos_z
                double_buffer[6*ii+3] = parts[i].core_part_space.vel_x
                double_buffer[6*ii+4] = parts[i].core_part_space.vel_y
                double_buffer[6*ii+5] = parts[i].core_part_space.vel_z
                labs[ii] = parts[i].lab
                ii = ii + 1
            end
        end
    elseif config[0].keytrj == 2 then
        var ii = 0
        for i in parts.ispace do
            if neighbour_init.check_valid(parts[i].neighbour_part_space) then
                double_buffer[9*ii] = parts[i].core_part_space.pos_x
                double_buffer[9*ii+1] = parts[i].core_part_space.pos_y
                double_buffer[9*ii+2] = parts[i].core_part_space.pos_z
                double_buffer[9*ii+3] = parts[i].core_part_space.vel_x
                double_buffer[9*ii+4] = parts[i].core_part_space.vel_y
                double_buffer[9*ii+5] = parts[i].core_part_space.vel_z
                double_buffer[9*ii+6] = parts[i].fxx
                double_buffer[9*ii+7] = parts[i].fyy
                double_buffer[9*ii+8] = parts[i].fzz
                labs[ii] = parts[i].lab
                ii = ii + 1
            end
        end
    end

   --Write the labs
    c_stdio.fwrite(labs, [terralib.sizeof(int32)], config[0].nsyst, HISTORY)
    --
    c_stdio.fwrite(double_buffer, [terralib.sizeof(double)],  3*(config[0].keytrj + 1) * config[0].nsyst, HISTORY)
    regentlib.c.free(labs)
    regentlib.c.free(double_buffer)
    c_stdio.fclose(HISTORY)
end

task write_export( config : region(ispace(int1d), config_type), parts : region(ispace(int1d), part), time : double) where 
                                                                                reads(config),
                                                                                        reads(parts.core_part_space,
                                                                                        parts.lab, parts.fxx, parts.fyy, parts.fzz,
                                                                                        parts.ltm, parts.ltp, parts.neighbour_part_space),
                                                                                writes(config.numframe, config.filesize)
do
    --Collect the variables
    var int_buf : &int32 = [&int32](regentlib.c.malloc( [terralib.sizeof( int32 )] * 3 * config[0].nsyst ))
    var double_buf : &double = [&double](regentlib.c.malloc( [terralib.sizeof( double)] * 9 * config[0].nsyst))
    regentlib.assert(not isnull(int_buf), "Failed to allocate double_buffer")
    regentlib.assert(not isnull(double_buf), "Failed to allocate double_buffer")
    var i : int32
    var ii : int32 = 0
    for i in parts.ispace do
        if neighbour_init.check_valid(parts[i].neighbour_part_space) then
            int_buf[ii*3] = parts[i].lab
            int_buf[ii*3+1] = parts[i].ltp
            int_buf[ii*3+2] = parts[i].ltm
            double_buf[9*ii] = parts[i].core_part_space.pos_x
            double_buf[9*ii+1] = parts[i].core_part_space.pos_y
            double_buf[9*ii+2] = parts[i].core_part_space.pos_z
            double_buf[9*ii+3] = parts[i].core_part_space.vel_x
            double_buf[9*ii+4] = parts[i].core_part_space.vel_y
            double_buf[9*ii+5] = parts[i].core_part_space.vel_z
            double_buf[9*ii+6] = parts[i].fxx
            double_buf[9*ii+7] = parts[i].fyy
            double_buf[9*ii+8] = parts[i].fzz
            ii = ii + 1
        end
    end
    
    var export = c_stdio.fopen('export', 'wb')
    var name : int8[81]
    c_string.strcpy(&(name[0]), config[0].text)
    c_stdio.fwrite(&(name[0]), [terralib.sizeof(int8)], 80, export)
    var nsyst = config[0].nsyst
    c_stdio.fwrite(&nsyst, [terralib.sizeof(int32)], 1, export)
    var nusyst = config[0].nusyst
    c_stdio.fwrite(&nusyst, [terralib.sizeof(int32)], 1, export)
    var temp = config[0].temp
    c_stdio.fwrite(&temp, [terralib.sizeof(double)], 1, export)
    var tstep = config[0].tstep
    c_stdio.fwrite(&tstep, [terralib.sizeof(double)], 1, export)
    var dimx = config[0].space.dim_x
    c_stdio.fwrite(&dimx, [terralib.sizeof(double)], 1, export)
    var dimy = config[0].space.dim_y
    c_stdio.fwrite(&dimy, [terralib.sizeof(double)], 1, export)
    var dimz = config[0].space.dim_z
    c_stdio.fwrite(&dimz, [terralib.sizeof(double)], 1, export)
    var shrdx = config[0].shrdx
    c_stdio.fwrite(&shrdx, [terralib.sizeof(double)], 1, export)
    var shrdy = config[0].shrdy
    c_stdio.fwrite(&shrdy, [terralib.sizeof(double)], 1, export)
    var shrdz = config[0].shrdz
    c_stdio.fwrite(&shrdz, [terralib.sizeof(double)], 1, export)

    c_stdio.fwrite(int_buf, [terralib.sizeof(int32)], 3 * config[0].nsyst, export)
    c_stdio.fwrite(double_buf, [terralib.sizeof(double)], 9 * config[0].nsyst, export)

    regentlib.c.free(int_buf)
    regentlib.c.free(double_buf)

    c_stdio.fclose(export)
 
end

__demand(__inline)
task write_revive( config : region(ispace(int1d), config_type) ) where reads(config) --.text,
                                                                            -- config.upx, config.upy, config.upz,
                                                                            -- config.fpx, config.fpy, config.fpz,
                                                                            -- config.nstep, config.nav, config.nstk,
                                                                            -- config.stpte, config.stppe, config.stpee,
                                                                            -- config.stpse, config.stpbe, config.stpae,
                                                                            -- config.stpde, config.stpvir, config.stptke,
                                                                            -- config.stpprs, config.stpvlm, config.stpzts,
                                                                            -- config.stpttp, config.stptpx, config.stptpy,
                                                                            -- config.stptpz, config.stpbdl, config.stpang, 
                                                                            -- config.stpdhd, config.rav, config.ave,
                                                                            -- config.flc, config.zum, config.nstk,
                                                                            -- config.stkpe, config.stkee, config.stkse,
                                                                            -- config.stkde, config.stkae, config.stkbe,
                                                                            -- config.stkvir, config.stkvlm, config.stkzts,
                                                                            -- config.stktkex, config.stktkey, config.stktkez,
                                                                            -- config.mt)
do
    var REVIVE = c_stdio.fopen('REVIVE', 'wb')
--     write simulation name, barostat properties and statistical accumulators
    var text : int8[81]
    c_string.strcpy(&(text[0]), config[0].text)
    c_stdio.fwrite(&(text[0]), [terralib.sizeof(int8)], 80, REVIVE)

    --Barostat properties
    var upx = config[0].upx
    var upy = config[0].upy
    var upz = config[0].upz
    var fpx = config[0].fpx
    var fpy = config[0].fpy
    var fpz = config[0].fpz
    c_stdio.fwrite(&upx, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&upy, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&upz, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&fpx, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&fpy, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&fpz, [terralib.sizeof(double)], 1, REVIVE)

    --timestep info
    var nstep = config[0].nstep
    var nav = config[0].nav
    var nstk = config[0].nstk
    c_stdio.fwrite(&nstep, [terralib.sizeof(int32)], 1, REVIVE)
    c_stdio.fwrite(&nav, [terralib.sizeof(int32)], 1, REVIVE)
    c_stdio.fwrite(&nstk, [terralib.sizeof(int32)], 1, REVIVE)

    --Statistics information
    var stpte = config[0].stpte
    var stppe = config[0].stppe
    var stpee = config[0].stpee
    var stpse = config[0].stpse
    var stpbe = config[0].stpbe
    var stpae = config[0].stpae
    var stpde = config[0].stpde
    c_stdio.fwrite(&stpte, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&stppe, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&stpee, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&stpse, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&stpbe, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&stpae, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&stpde, [terralib.sizeof(double)], 1, REVIVE)
    var stpvir = config[0].stpvir
    var stptke = config[0].stptke
    var stpprs = config[0].stpprs
    var stpvlm = config[0].stpvlm
    var stpzts = config[0].stpzts
    var stpttp = config[0].stpttp
    c_stdio.fwrite(&stpvir, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&stptke, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&stpprs, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&stpvlm, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&stpzts, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&stpttp, [terralib.sizeof(double)], 1, REVIVE)
    var stptpx = config[0].stptpx
    var stptpy = config[0].stptpy
    var stptpz = config[0].stptpz
    c_stdio.fwrite(&stptpx, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&stptpy, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&stptpz, [terralib.sizeof(double)], 1, REVIVE)
    var stpbdl = config[0].stpbdl
    var stpang = config[0].stpang
    var stpdhd = config[0].stpdhd
    c_stdio.fwrite(&stpbdl, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&stpang, [terralib.sizeof(double)], 1, REVIVE)
    c_stdio.fwrite(&stpdhd, [terralib.sizeof(double)], 1, REVIVE)
    var rav : double[CONST_stksize]
    for i=0, CONST_stksize do
        rav[i] = config[0].rav[i]
    end
    c_stdio.fwrite(&(rav[0]), [terralib.sizeof(double)], CONST_stksize, REVIVE)
    var ave : double[CONST_statsize]
    var flc : double[CONST_statsize]
    for i = 0, CONST_statsize do
        ave[i] = config[0].ave[i]
        flc[i] = config[0].flc[i]
    end
    c_stdio.fwrite(&(ave[0]), [terralib.sizeof(double)], CONST_statsize, REVIVE)
    c_stdio.fwrite(&(flc[0]), [terralib.sizeof(double)], CONST_statsize, REVIVE)
    var zum : double[13]
    for i = 0, 13 do
        zum[i] = config[0].zum[i]
    end
    c_stdio.fwrite(&(zum[0]), [terralib.sizeof(double)], 13, REVIVE)
    var stkpe : double[CONST_maxstk]
    var stkee : double[CONST_maxstk]
    var stkse : double[CONST_maxstk]
    var stkde : double[CONST_maxstk]
    var stkae : double[CONST_maxstk]
    var stkbe : double[CONST_maxstk]
    for i = 0, config[0].nstk do
        stkpe[i] = config[0].stkpe[i]
        stkee[i] = config[0].stkee[i]
        stkse[i] = config[0].stkse[i]
        stkde[i] = config[0].stkde[i]
        stkae[i] = config[0].stkae[i]
        stkbe[i] = config[0].stkbe[i]
    end
    c_stdio.fwrite(&(stkpe[0]), [terralib.sizeof(double)], config[0].nstk, REVIVE)
    c_stdio.fwrite(&(stkee[0]), [terralib.sizeof(double)], config[0].nstk, REVIVE)
    c_stdio.fwrite(&(stkse[0]), [terralib.sizeof(double)], config[0].nstk, REVIVE)
    c_stdio.fwrite(&(stkde[0]), [terralib.sizeof(double)], config[0].nstk, REVIVE)
    c_stdio.fwrite(&(stkae[0]), [terralib.sizeof(double)], config[0].nstk, REVIVE)
    c_stdio.fwrite(&(stkbe[0]), [terralib.sizeof(double)], config[0].nstk, REVIVE)
    var stkvir : double[CONST_maxstk]
    var stkvlm : double[CONST_maxstk]
    var stkzts : double[CONST_maxstk]
    var stktkex: double[CONST_maxstk]
    var stktkey: double[CONST_maxstk]
    var stktkez: double[CONST_maxstk]
    for i = 0, config[0].nstk do
        stkvir[i] = config[0].stkvir[i]
        stkvlm[i] = config[0].stkvlm[i]
        stkzts[i] = config[0].stkzts[i]
        stktkex[i] = config[0].stktkex[i]
        stktkey[i] = config[0].stktkey[i]
        stktkez[i] = config[0].stktkez[i]
    end
    c_stdio.fwrite(&(stkvir[0]), [terralib.sizeof(double)], config[0].nstk, REVIVE)
    c_stdio.fwrite(&(stkvlm[0]), [terralib.sizeof(double)], config[0].nstk, REVIVE)
    c_stdio.fwrite(&(stkzts[0]), [terralib.sizeof(double)], config[0].nstk, REVIVE)
    c_stdio.fwrite(&(stktkex[0]), [terralib.sizeof(double)], config[0].nstk, REVIVE)
    c_stdio.fwrite(&(stktkey[0]), [terralib.sizeof(double)], config[0].nstk, REVIVE)
    c_stdio.fwrite(&(stktkez[0]), [terralib.sizeof(double)], config[0].nstk, REVIVE)

    var nodes : int32 = 1
    var threads : int32 = 1
    --TODO: NYI: Checking threads
    c_stdio.fwrite(&nodes, [terralib.sizeof(int32)], 1, REVIVE)
    c_stdio.fwrite(&threads, [terralib.sizeof(int32)], 1, REVIVE)
   
    --TODO : Random state...?
    var mt : int32[CONST_randsize]
    for i = 0, CONST_randsize do
        mt[i] = config[0].mt[i]
    end
    c_stdio.fwrite(&(mt[0]), [terralib.sizeof(int32)], CONST_randsize, REVIVE)
 
    c_stdio.fclose(REVIVE)
end


local sqrtf = regentlib.sqrt(double)

--TODO: NYI write_config
task write_config()
    regentlib.assert(false, "write_config not yet implemented.")
end

task write_output_result( config : region(ispace(int1d), config_type), parts : region(ispace(int1d), part) ) where
                                                                        reads(config), --.text,
--                                                                             config.upx, config.upy, config.upz,
--                                                                             config.fpx, config.fpy, config.fpz,
--                                                                             config.nstep, config.nav, config.nstk,
--                                                                             config.stpte, config.stppe, config.stpee,
--                                                                             config.stpse, config.stpbe, config.stpae,
--                                                                             config.stpde, config.stpvir, config.stptke,
--                                                                             config.stpprs, config.stpvlm, config.stpzts,
--                                                                             config.stpttp, config.stptpx, config.stptpy,
--                                                                             config.stptpz, config.stpbdl, config.stpang, 
--                                                                             config.stpdhd, config.rav, config.ave,
--                                                                             config.flc, config.zum, config.nstk,
--                                                                             config.stkpe, config.stkee, config.stkse,
--                                                                             config.stkde, config.stkae, config.stkbe,
--                                                                             config.stkvir, config.stkvlm, config.stkzts,
--                                                                             config.stktkex, config.stktkey, config.stktkez,
--                                                                             config.mt), reads(config.nsyst),
                                                                                        reads(parts.core_part_space,
                                                                                        parts.lab, parts.atmnam),
                                                                        writes(config.flc)
do
  --Line 3044
    write_revive(config)
    var OUTPUT = c_stdio.fopen("OUTPUT", "a")
    var i : int32
    
    var buffer : int8[25]
    c_stdio.snprintf(&(buffer[0]), 25, "%i", config[0].nstep)
    var buffer2 : int8[25]
    c_stdio.snprintf(&(buffer2[0]), 25, "%i", config[0].nav)

    c_stdio.fprintf(OUTPUT, "\n run closing at step %10s final averages and fluctuations over %10s steps", buffer, buffer2)
    for i = 0, CONST_statsize do
        config[0].flc[i] = sqrtf(config[0].flc[i])
    end

    --Assumed for now, config[0].outsel is 0
    --if config[0].outsel == 0 then
        c_stdio.fprintf(OUTPUT, "\n ")
        for i=0, 95 do
            c_stdio.fprintf(OUTPUT, "-")
        end
        c_stdio.fprintf(OUTPUT, "\n       ")
        c_stdio.fprintf(OUTPUT, "en-total      pe-total     vir-total      ke-total      pressure   temperature\n ")
        for i=0, 95 do
            c_stdio.fprintf(OUTPUT, "-")
        end
        c_stdio.fprintf(OUTPUT, "\n") 
        c_stdio.fprintf(OUTPUT, "  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", config[0].ave[0], config[0].ave[1],
                                config[0].ave[7], config[0].ave[8], config[0].ave[9], config[0].ave[12] )
        c_stdio.fprintf(OUTPUT, "  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", config[0].flc[0], config[0].flc[1],
                                config[0].flc[7], config[0].flc[8], config[0].flc[9], config[0].flc[12] )
        for i=0, 95 do
            c_stdio.fprintf(OUTPUT, "-")
        end
        c_stdio.fprintf(OUTPUT, "\n") 
    --elseif...
    --end
    --write out average stress tensors (conservative, dissipative, random and kinetic)
    c_stdio.fprintf(OUTPUT, "\n\n         average conservative stress tensor                       r.m.s. fluctuations\n\n")
    c_stdio.fprintf(OUTPUT, "  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", config[0].ave[16], config[0].ave[20],
                            config[0].ave[24], config[0].flc[16], config[0].flc[20], config[0].flc[24] )
    c_stdio.fprintf(OUTPUT, "  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", config[0].ave[28], config[0].ave[32],
                            config[0].ave[36], config[0].flc[28], config[0].flc[32], config[0].flc[36] )
    c_stdio.fprintf(OUTPUT, "  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", config[0].ave[40], config[0].ave[44],
                            config[0].ave[48], config[0].flc[40], config[0].flc[44], config[0].flc[48] )

    c_stdio.fprintf(OUTPUT, "\n\n          average dissipative stress tensor                       r.m.s. fluctuations\n\n")
    c_stdio.fprintf(OUTPUT, "  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", config[0].ave[17], config[0].ave[21],
                            config[0].ave[25], config[0].flc[17], config[0].flc[21], config[0].flc[25] )
    c_stdio.fprintf(OUTPUT, "  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", config[0].ave[29], config[0].ave[33],
                            config[0].ave[37], config[0].flc[29], config[0].flc[33], config[0].flc[37] )
    c_stdio.fprintf(OUTPUT, "  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", config[0].ave[41], config[0].ave[45],
                            config[0].ave[49], config[0].flc[41], config[0].flc[45], config[0].flc[49] )

    c_stdio.fprintf(OUTPUT, "\n\n               average random stress tensor                       r.m.s. fluctuations\n\n")
    c_stdio.fprintf(OUTPUT, "  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", config[0].ave[18], config[0].ave[22],
                            config[0].ave[26], config[0].flc[18], config[0].flc[22], config[0].flc[26] )
    c_stdio.fprintf(OUTPUT, "  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", config[0].ave[30], config[0].ave[34],
                            config[0].ave[38], config[0].flc[30], config[0].flc[34], config[0].flc[38] )
    c_stdio.fprintf(OUTPUT, "  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", config[0].ave[42], config[0].ave[46],
                            config[0].ave[50], config[0].flc[42], config[0].flc[46], config[0].flc[50] )

    c_stdio.fprintf(OUTPUT, "\n\n              average kinetic stress tensor                       r.m.s. fluctuations\n\n")
    c_stdio.fprintf(OUTPUT, "  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", config[0].ave[19], config[0].ave[23],
                            config[0].ave[27], config[0].flc[19], config[0].flc[23], config[0].flc[27] )
    c_stdio.fprintf(OUTPUT, "  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", config[0].ave[31], config[0].ave[35],
                            config[0].ave[39], config[0].flc[31], config[0].flc[35], config[0].flc[39] )
    c_stdio.fprintf(OUTPUT, "  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", config[0].ave[43], config[0].ave[47],
                            config[0].ave[51], config[0].flc[43], config[0].flc[47], config[0].flc[51] )

    --Not yet recording runtime here
    c_stdio.fprintf(OUTPUT, "\n average cpu time (forces) = %10.5f (s)\n average cpu time (cycle)  = %10.5f (s)\n", 0.0, 0.0)

    --TODO: NYI final part/max interaction counters
    c_stdio.fprintf(OUTPUT, "\n final no. buffer particles   = %14s\n final max. interactions/node = %14s\n", "0", "0")
    
    if(config[0].nsyst > 0) then
        c_stdio.fprintf(OUTPUT, "\n final particle positions and velocities\n\n")
        var k = regentlib.fmax(1, config[0].nsyst/20)
        for i=0, config[0].nsyst, k do --FIXME: Validity check
            --if bond then
            --else
                c_stdio.snprintf(&(buffer[0]), 25, "%i", parts[int1d(i)].lab)
                c_stdio.fprintf(OUTPUT, "%10s %8s %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",
                                        buffer, parts[int1d(i)].atmnam, parts[int1d(i)].core_part_space.pos_x,   
                                        parts[int1d(i)].core_part_space.pos_y, parts[int1d(i)].core_part_space.pos_z,
                                        parts[int1d(i)].core_part_space.vel_x, parts[int1d(i)].core_part_space.vel_y,
                                        parts[int1d(i)].core_part_space.vel_z  )
            --end
        end
    end

    c_stdio.fclose(OUTPUT)
end

task gather_write_data( lexport : bool, lhistory : bool, time : double, 
                        config : region(ispace(int1d), config_type), parts : region(ispace(int1d), part)) where 
                                                                                reads(config),
                                                                                        reads(parts.core_part_space,
                                                                                        parts.lab, parts.fxx, parts.fyy, parts.fzz,
                                                                                        parts.ltm, parts.ltp, parts.neighbour_part_space),
                                                                                writes(config.numframe, config.filesize)

do

    if ( lhistory) then
        write_history(config ,parts , time)
    end

    if( lexport) then
        write_export(config, parts, time)
    end

end

task write_output_header(config : region(ispace(int1d), config_type)) where reads(config.l_scr)
do
    regentlib.assert(config[0].l_scr == false, "Not yet supporting output to stdio")
    
    var OUTPUT = c_stdio.fopen('OUTPUT', 'w') --Overwrites any existing file
    c_stdio.fprintf(OUTPUT, " ***********************  HartreeParticleDSL implementation of  \n")
    c_stdio.fprintf(OUTPUT, "                          ukri stfc daresbury laboratory dissipative\n")
    c_stdio.fprintf(OUTPUT, "                          particle dynamics program DL_MESO_DPD\n")
    c_stdio.fprintf(OUTPUT, "                          author: a.b.g. chalk\n")
    c_stdio.fprintf(OUTPUT, "  D L _ M E S O _ D P D   original authors of DL_MESO_DPD:\n")
    c_stdio.fprintf(OUTPUT, "                                        w. smith & m. a. seaton\n")
    c_stdio.fprintf(OUTPUT, "                                                                   \n")
    c_stdio.fprintf(OUTPUT, "                          Based upon dl_meso version:\n")
    c_stdio.fprintf(OUTPUT, "                          dl_meso version 2.8 rev 00, october 2020\n")
    c_stdio.fprintf(OUTPUT, " ***********************  copyright (c) 2021 UKRI STFC Hartree Centre\n\n")

    c_stdio.fprintf(OUTPUT, " **************************************************************************\n")
    c_stdio.fprintf(OUTPUT, " *********************       RUNNING WITH REGENT        *******************\n")
    c_stdio.fprintf(OUTPUT, " *********************         (LITTLE ENDIAN)          *******************\n")
    c_stdio.fprintf(OUTPUT, " **************************************************************************\n")
   
    c_stdio.fclose(OUTPUT) 
end

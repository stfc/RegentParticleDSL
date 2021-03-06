-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

local privilege_lists = {}

--Compute the privileges used for a symmetric pair task.
-- Inputs:
-- parts1: Symbol used for the first particle region in the task.
-- parts2: Symbol used for the second particle region in the task.
-- read1/read2: Output from the compute_privileges code for the fields read.
-- write1/write2: Output from the compute_privileges code for the field written to.
-- reduc1/reduc2: Output from the compute_privileges code for the fields reduced.
function privilege_lists.get_privileges_symmetric_pair_task( parts1, parts2, config, read1, read2, read3, write1, write2, reduc1, reduc2, reduc3 )
    local read1_privs = terralib.newlist()
    local read2_privs = terralib.newlist()
    local read3_privs = terralib.newlist()
    local write1_privs = terralib.newlist()
    local write2_privs = terralib.newlist()
    local reduc1_privs = terralib.newlist()
    local reduc2_privs = terralib.newlist()
    local reduc3_privs = terralib.newlist()
    local update_neighbours = false

    for _, v in pairs(read1) do
        read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
    end
    for _, v in pairs(read2) do
        read2_privs:insert( regentlib.privilege(regentlib.reads, parts2, string_to_field_path.get_field_path(v)))
    end
    for _, v in pairs(read3) do
        read3_privs:insert( regentlib.privilege(regentlib.reads, config, string_to_field_path.get_field_path(v)))
    end
    for _, v in pairs(write1) do
        write1_privs:insert( regentlib.privilege(regentlib.writes, parts1, string_to_field_path.get_field_path(v)))
        if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" then
            update_neighbours = true
        end
    end
    for _, v in pairs(write2) do
        write2_privs:insert( regentlib.privilege(regentlib.writes, parts2, string_to_field_path.get_field_path(v)))
        if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" then
            update_neighbours = true
        end
    end
    for _, v in pairs(reduc1) do
        reduc1_privs:insert( regentlib.privilege(regentlib.reduces(v[2]), parts1, string_to_field_path.get_field_path(v[1])))
        if v[1] == "core_part_space.pos_x" or v[1] == "core_part_space.pos_y" or v[1] == "core_part_space.pos_z" then
            update_neighbours = true
        end
    end
    for _, v in pairs(reduc2) do
        reduc2_privs:insert( regentlib.privilege(regentlib.reduces(v[2]), parts2, string_to_field_path.get_field_path(v[1])))
        if v[1] == "core_part_space.pos_x" or v[1] == "core_part_space.pos_y" or v[1] == "core_part_space.pos_z" then
            update_neighbours = true
        end
    end
    for _, v in pairs(reduc3) do
        reduc3_privs:insert( regentlib.privilege(regentlib.reduces(v[2]), config, string_to_field_path.get_field_path(v[1])))
    end


    return update_neighbours, read1_privs, read2_privs, read3_privs, write1_privs, write2_privs, reduc1_privs, reduc2_privs, reduc3_privs
end

function privilege_lists.get_privileges_asymmetric_pair_task( parts1, parts2, config, read1, read2,read3, write1, reduc1, reduc3 )

    local read1_privs = terralib.newlist()
    local read2_privs = terralib.newlist()
    local read3_privs = terralib.newlist()
    local write1_privs = terralib.newlist()
    local reduc1_privs = terralib.newlist()
    local reduc3_privs = terralib.newlist()
    local update_neighbours = false
    
    for _, v in pairs(read1) do
        read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
    end
    for _, v in pairs(read2) do
        read2_privs:insert( regentlib.privilege(regentlib.reads, parts2, string_to_field_path.get_field_path(v)))
    end  
    for _, v in pairs(read3) do
        read3_privs:insert( regentlib.privilege(regentlib.reads, config, string_to_field_path.get_field_path(v)))
    end
    for _, v in pairs(write1) do
        write1_privs:insert( regentlib.privilege(regentlib.writes, parts1, string_to_field_path.get_field_path(v)))
        if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" then
            update_neighbours = true
        end
    end
    for _, v in pairs(reduc1) do
        reduc1_privs:insert( regentlib.privilege(regentlib.reduces(v[2]), parts1, string_to_field_path.get_field_path(v[1])))
        if v[1] == "core_part_space.pos_x" or v[1] == "core_part_space.pos_y" or v[1] == "core_part_space.pos_z" then
            update_neighbours = true
        end
    end
    for _, v in pairs(reduc3) do
        reduc3_privs:insert( regentlib.privilege(regentlib.reduces(v[2]), config, string_to_field_path.get_field_path(v[1])))
    end

    return update_neighbours, read1_privs, read2_privs, read3_privs, write1_privs, reduc1_privs, reduc3_privs
end

function privilege_lists.get_privileges_self_task( parts1, config, read1, read2, read3, write1, write2, reduc1, reduc2, reduc3 )

    local read1_privs = terralib.newlist()
    local read2_privs = terralib.newlist()
    local write1_privs = terralib.newlist()
    local reduc1_privs = terralib.newlist()
    local reduc2_privs = terralib.newlist()
    local update_neighbours = false
    
    for _, v in pairs(read1) do
        read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
    end
    for _, v in pairs(read2) do
        read1_privs:insert( regentlib.privilege(regentlib.reads, parts1, string_to_field_path.get_field_path(v)))
    end
    for _, v in pairs(read3) do
        read2_privs:insert(regentlib.privilege(regentlib.reads, config, string_to_field_path.get_field_path(v)))
    end
    for _, v in pairs(write1) do
        write1_privs:insert( regentlib.privilege(regentlib.writes, parts1, string_to_field_path.get_field_path(v)))
        if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" then
            update_neighbours = true
        end
    end
    for _, v in pairs(write2) do
        write1_privs:insert( regentlib.privilege(regentlib.writes, parts1, string_to_field_path.get_field_path(v)))
        if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" then
            update_neighbours = true
        end
    end
    for _, v in pairs(reduc1) do
        reduc1_privs:insert( regentlib.privilege(regentlib.reduces(v[2]), parts1, string_to_field_path.get_field_path(v[1])))
        if v[1] == "core_part_space.pos_x" or v[1] == "core_part_space.pos_y" or v[1] == "core_part_space.pos_z" then
            update_neighbours = true
        end
    end
    for _, v in pairs(reduc2) do
        reduc1_privs:insert( regentlib.privilege(regentlib.reduces(v[2]), parts1, string_to_field_path.get_field_path(v[1])))
        if v[1] == "core_part_space.pos_x" or v[1] == "core_part_space.pos_y" or v[1] == "core_part_space.pos_z" then
            update_neighbours = true
        end
    end
    for _, v in pairs(reduc3) do
        reduc2_privs:insert( regentlib.privilege(regentlib.reduces(v[2]), config, string_to_field_path.get_field_path(v[1])))
    end
    return update_neighbours, read1_privs, read2_privs, write1_privs, reduc1_privs, reduc2_privs
end

function privilege_lists.get_privileges_per_part( part1, read1, write1, reduc1, config, configread, configreduc)

    local read1_privs = terralib.newlist()
    local write1_privs = terralib.newlist()
    local reduc1_privs = terralib.newlist()
    local readconf_privs = terralib.newlist()
    local reducconf_privs = terralib.newlist()
    local update_neighbours = false
    for _, v in pairs(read1) do
        read1_privs:insert( regentlib.privilege(regentlib.reads, part1, string_to_field_path.get_field_path(v)))
    end
    for _, v in pairs(write1) do
        write1_privs:insert( regentlib.privilege(regentlib.writes, part1, string_to_field_path.get_field_path(v)))
        if v == "core_part_space.pos_x" or v == "core_part_space.pos_y" or v == "core_part_space.pos_z" then
            update_neighbours = true
        end
    end
    for _, v in pairs(reduc1) do
        reduc1_privs:insert( regentlib.privilege(regentlib.reduces(v[2]), part1, string_to_field_path.get_field_path(v[1])))
        if v[1] == "core_part_space.pos_x" or v[1] == "core_part_space.pos_y" or v[1] == "core_part_space.pos_z" then
            update_neighbours = true
        end
    end
    for _,v in pairs(configread) do
        readconf_privs:insert( regentlib.privilege(regentlib.reads, config, string_to_field_path.get_field_path(v)))
    end
    for _, v in pairs(configreduc) do
        reducconf_privs:insert( regentlib.privilege(regentlib.reduces(v[2]), config, string_to_field_path.get_field_path(v[1])))
    end

    return update_neighbours, read1_privs, write1_privs, reduc1_privs, readconf_privs, reducconf_privs
end

return privilege_lists


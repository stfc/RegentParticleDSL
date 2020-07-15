import "regent"

require("example_types/example_interaction")

local c = regentlib.c
local format = require("std/format")

function init_two_particles()

local init_string = terralib.newlist()

init_string:insert(rquote
    var particle_one = region(ispace(int1d, 1), part)
    var particle_two = region(ispace(int1d, 1), part)
end)


return init_string
end

function make_pairwise_interaction(kernel_name, function_name)

local parts1 = regentlib.newsymbol(region(ispace(int1d),part), "parts1")
local parts2 = regentlib.newsymbol(region(ispace(int1d),part), "parts2")
local part1 = regentlib.newsymbol("part1")
local part2 = regentlib.newsymbol("part2")

local interaction = kernel_name(part1, part2)

local task_name = regentlib.newsymbol("pairwise_"..function_name)
local task pairwise_task(parts1 : region(ispace(int1d),part), parts2 : region(ispace(int1d),part)) where reads(parts1, parts2), writes(parts1, parts2) do
   for [part1] in [parts1] do
     for [part2] in [parts2] do
       [parts1][part1].extra_variable_2 = 50
       [interaction]
     end
   end
end

return pairwise_task
end

local pairwise_task = make_pairwise_interaction(kernel_one, "one")
local pairwise_task_two = make_pairwise_interaction(kernel_two, "two")

task main()
    var particle_one = region(ispace(int1d, 1), part)
    var particle_two = region(ispace(int1d, 1), part)
    fill(particle_one.{neighbour_part_space.cell_id}, int3d({0,0,0}))
    fill(particle_one.{core_part_space.pos_x}, 0) 
    fill(particle_one.{core_part_space.pos_y}, 0) 
    fill(particle_one.{core_part_space.pos_z}, 0) 
    fill(particle_one.{core_part_space.mass}, 0) 
    fill(particle_one.{core_part_space.vel_x}, 0) 
    fill(particle_one.{core_part_space.vel_y}, 0) 
    fill(particle_one.{core_part_space.vel_z}, 0) 
    fill(particle_one.{core_part_space.cutoff}, 0) 
    fill(particle_one.{core_part_space.id}, 0) 
    fill(particle_one.{extra_variable_1}, 0) 
    fill(particle_one.{extra_variable_2}, 0) 
    fill(particle_one.{extra_variable_3}, 0) 
    fill(particle_two.{neighbour_part_space.cell_id}, int3d({0,0,0}))
    fill(particle_two.{core_part_space.pos_x}, 0) 
    fill(particle_two.{core_part_space.pos_y}, 0) 
    fill(particle_two.{core_part_space.pos_z}, 0) 
    fill(particle_two.{core_part_space.mass}, 0) 
    fill(particle_two.{core_part_space.vel_x}, 0) 
    fill(particle_two.{core_part_space.vel_y}, 0) 
    fill(particle_two.{core_part_space.vel_z}, 0) 
    fill(particle_two.{core_part_space.cutoff}, 0) 
    fill(particle_two.{core_part_space.id}, 0) 
    fill(particle_two.{extra_variable_1}, 0) 
    fill(particle_two.{extra_variable_2}, 0) 
    fill(particle_two.{extra_variable_3}, 0)
    pairwise_task(particle_one, particle_two)
    pairwise_task_two(particle_one, particle_two)
    c.legion_runtime_issue_execution_fence(__runtime(), __context())

    format.println("{}",particle_one[0].extra_variable_1)
    format.println("{}", particle_one[0].extra_variable_2)
    format.println("{}",particle_two[0].extra_variable_1) 
end

regentlib.start(main)

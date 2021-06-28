import "regent"
require("src/RegentParticleDSL")
set_dimensionality(3)
set_periodicity(true)
disable_high_performance()
enable_timing()
setup_part()
local format = require("std/format")

fspace part{
    core_part_space : core_part,
    neighbour_part_space : neighbour_part,
    a : int32[4]
}
setup_dsl()
require("src/io_modules/empty_io/import_empty_io")

function mover(part1, config)
    local kernel = rquote
        part1.core_part_space.pos_x = part1.core_part_space.pos_x +  1.0
    end
    return kernel
end


local c_stdio = terralib.includec("stdio.h")

task main()

[initialisation_function(variables, 2, 2.5, 2.5, 2.5)];
variables.particle_array[0].a[0] = 34;
variables.particle_array[0].a[1] = 512;
variables.particle_array[0].a[2] = 256;
variables.particle_array[0].core_part_space.cutoff = 0.5
variables.particle_array[0].core_part_space.pos_x = 0.5
variables.particle_array[0].core_part_space.pos_y = 0.5
variables.particle_array[0].core_part_space.pos_z = 0.5

variables.particle_array[1].a[0] = 34;
variables.particle_array[1].a[1] = 512;
variables.particle_array[1].a[2] = 256;
variables.particle_array[1].core_part_space.cutoff = 0.5
variables.particle_array[1].core_part_space.pos_x = 0.5
variables.particle_array[1].core_part_space.pos_y = 1.0
variables.particle_array[1].core_part_space.pos_z = 0.5


[neighbour_init.initialise(variables)];
[neighbour_init.update_cells(variables)];

for part in neighbour_init.padded_particle_array.ispace do
   if not neighbour_init.check_valid(neighbour_init.padded_particle_array[part].neighbour_part_space) then
        neighbour_init.padded_particle_array[part].a[0] = 123;
        neighbour_init.padded_particle_array[part].a[1] = 123;
        neighbour_init.padded_particle_array[part].a[2] = 123;
        neighbour_init.padded_particle_array[part].core_part_space.pos_z = 0.0
    end
end


[invoke(variables.config, {mover, PER_PART}, BARRIER)];

--for part in neighbour_init.padded_particle_array.ispace do
--   if not neighbour_init.check_valid(neighbour_init.padded_particle_array[part].neighbour_part_space) then
--        neighbour_init.padded_particle_array[part].a[0] = 123;
--        neighbour_init.padded_particle_array[part].a[1] = 123;
--        neighbour_init.padded_particle_array[part].a[2] = 123;
--    end
--end

for part in neighbour_init.padded_particle_array do
    if neighbour_init.check_valid(neighbour_init.padded_particle_array[part].neighbour_part_space) then
        regentlib.assert(neighbour_init.padded_particle_array[part].a[0] == 34, "Failed to find value 34 in a[0]")
        regentlib.assert(neighbour_init.padded_particle_array[part].a[1] == 512,"Failed to find value 512 in a[1]")
        regentlib.assert(neighbour_init.padded_particle_array[part].a[2] == 256, "Failed to find value 256 in a[2]")
    else
        regentlib.assert(neighbour_init.padded_particle_array[part].a[0] == 123 or 
                        neighbour_init.padded_particle_array[part].core_part_space.pos_z== 0.5, "Value wrong in invalid particle")
        regentlib.assert(neighbour_init.padded_particle_array[part].a[1] == 123 or
                        neighbour_init.padded_particle_array[part].core_part_space.pos_z== 0.5, "Value wrong in invalid particle")
        regentlib.assert(neighbour_init.padded_particle_array[part].a[2] == 123 or
                         neighbour_init.padded_particle_array[part].core_part_space.pos_z== 0.5, "Value wrong in invalid particle")
    end
end

for part in neighbour_init.padded_particle_array.ispace do
    if neighbour_init.check_valid(neighbour_init.padded_particle_array[part].neighbour_part_space) then
        format.println("Setting {}", part)
        neighbour_init.padded_particle_array[part].a[0] = 13
        neighbour_init.padded_particle_array[part].a[1] = 13
        neighbour_init.padded_particle_array[part].a[2] = 13
    end
end


for part in neighbour_init.padded_particle_array.ispace do
    if (not neighbour_init.check_valid(neighbour_init.padded_particle_array[part].neighbour_part_space)) and neighbour_init.padded_particle_array[part].core_part_space.pos_z== 0.5 then
        format.println("Checking {}",part)
        regentlib.assert(neighbour_init.padded_particle_array[part].a[0] == 34, "Failed at the side-effect check")
        regentlib.assert(neighbour_init.padded_particle_array[part].a[1] == 512, "Failed at the side-effect check")
        regentlib.assert(neighbour_init.padded_particle_array[part].a[2] == 256, "Failed at the side-effect check")
    end
end

end

run_DSL(main)

import "regent"

local format = require("std/format")

local M_PI = 3.141592653589793238462643383279502884
dx = global(double, 0.0025)
h = global(double, 0.0025 * 1.3)
h2 = global(double, h:get()*h:get())
h3 = global(double, h2:get()*h:get())
ad_7 = global(double, 3.0/(359.0*M_PI*h3:get()))
ad_7h = global(double, ad_7:get() / h:get())

require("src/RegentParticleDSL")
set_dimensionality(2)
set_periodicity(false)
setup_part()
local format = require("std/format")
require("src/interactions/ISPH/ISPH_part")
setup_dsl()
isph_module = require("src/io_modules/ISPH/2d_isph_module")

--variables = {}
--variables.config = regentlib.newsymbol("config")
--variables.particle_array = regentlib.newsymbol("particle_region")

require("src/interactions/ISPH/divergence")

task set_cutoff(particles : region(ispace(int1d), part)) where reads(particles), writes(particles) do

  for part in particles.ispace do
    particles[part].core_part_space.cutoff = 3.0 * [h:get()]
    particles[part].divergence = 0.0
    particles[part].core_part_space.id = part
  end
  particles[0].divergence = 0.0
end

task main()
[isph_module.initialisation_function("/home/aidan/isph_read/rec000000.txt", variables, 4.0, 3.0)];
set_cutoff(variables.particle_array);
[neighbour_init.initialise(variables)];
--[neighbour_init.update_cells(variables)];

[invoke(variables.config, {divergence,ASYMMETRIC_PAIRWISE}, NO_BARRIER)];

format.println("{} {}", neighbour_init.padded_particle_array[0].divergence, neighbour_init.padded_particle_array[357558].divergence)
format.println("{} {}", neighbour_init.padded_particle_array[0].interactions, neighbour_init.padded_particle_array[357558].interactions)
end

run_DSL(main)

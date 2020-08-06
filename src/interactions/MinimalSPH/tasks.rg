import "regent"

require("defaults")
--require("src/neighbour_search/cell_pair/import_cell_pair")
require("src/neighbour_search/cell_pair/neighbour_search")
require("src/neighbour_search/cell_pair/cell")
require("src/interactions/MinimalSPH/interactions")
require("src/interactions/MinimalSPH/timestep")

local density_task = create_asymmetric_pairwise_runner(nonsym_density_kernel)
local c = regentlib.c

--Some big number
local hmax = 12345678.0
local hmin = 0.0
--Rough estimate at the moment
local kernel_root = 0.4184291064739227294921875
local kernel_dimension = 3.0
local hydro_dimension_inv = 1.0 / kernel_dimension
local hydro_eta = 1.2348
local hydro_eps = 1e-4
local alpha = 0.8

local fabsd = regentlib.fabs(double)
local cbrtd = regentlib.cbrt(double)
local sqrt = regentlib.sqrt(double)
--local fmind = regentlib.fmin()
--local fmaxd = regentlib.fmax()
--3D case
terra cube_val(value : double)
  var x2 = value * value
  return x2 * value
end

terra fmind(val1 : double, val2 : double)
  return regentlib.fmin(val1, val2)
end

terra fmaxd(val1 : double, val2 : double)
  return regentlib.fmax(val1, val2)
end

local run_subset_task = run_asymmetric_subset_task(nonsym_density_kernel)

task pair_redo_density( parts_self : region(ispace(int1d), part), 
                        subset: partition(disjoint, parts_self, ispace(int1d)), 
                        parts_far : region(ispace(int1d), part),
                        space : region(ispace(int1d), space_config) )
   where reads(parts_self, parts_far, space), writes( parts_self ) do
   var no_redo = int1d(1)
   var redo = int1d(0)
   var half_box_x = 0.5 * space[0].dim_x
   var half_box_y = 0.5 * space[0].dim_y
   var half_box_z = 0.5 * space[0].dim_z
   for part1 in subset[redo].ispace do
     for part2 in parts_far.ispace do
       --Compute particle distance
         var dx = subset[redo][part1].core_part_space.pos_x - parts_far[part2].core_part_space.pos_x
         var dy = subset[redo][part1].core_part_space.pos_y - parts_far[part2].core_part_space.pos_y
         var dz = subset[redo][part1].core_part_space.pos_z - parts_far[part2].core_part_space.pos_z
         if (dx > half_box_x) then dx = dx - half_box_x end
         if (dy > half_box_y) then dy = dy - half_box_y end
         if (dz > half_box_z) then dz = dz - half_box_z end
         if (dx <-half_box_x) then dx = dx + half_box_x end
         if (dy <-half_box_y) then dy = dy + half_box_y end
         if (dz <-half_box_z) then dz = dz + half_box_z end
         var cutoff2 = subset[redo][part1].core_part_space.cutoff
         cutoff2 = cutoff2 * cutoff2
         var r2 = dx*dx + dy*dy + dz*dz
         if(r2 <= cutoff2) then
           [nonsym_density_kernel(rexpr subset[redo][part1] end, rexpr parts_far[part2] end, rexpr r2 end)]
         end
     end
   end
end


task self_redo_density( parts : region(ispace(int1d), part), subset : partition(disjoint, parts, ispace(int1d)), space : region(ispace(int1d), space_config) )
   where reads(parts, space), writes(parts) do
   var no_redo = int1d(1)
   var redo = int1d(0)
   var half_box_x = 0.5 * space[0].dim_x
   var half_box_y = 0.5 * space[0].dim_y
   var half_box_z = 0.5 * space[0].dim_z
   for part1 in subset[redo].ispace do
     for part2 in parts.ispace do
       --Compute particle distance
       if(int1d(part1) ~= int1d(part2)) then
         var dx = subset[redo][part1].core_part_space.pos_x - parts[part2].core_part_space.pos_x
         var dy = subset[redo][part1].core_part_space.pos_y - parts[part2].core_part_space.pos_y
         var dz = subset[redo][part1].core_part_space.pos_z - parts[part2].core_part_space.pos_z
         if (dx > half_box_x) then dx = dx - half_box_x end
         if (dy > half_box_y) then dy = dy - half_box_y end
         if (dz > half_box_z) then dz = dz - half_box_z end
         if (dx <-half_box_x) then dx = dx + half_box_x end
         if (dy <-half_box_y) then dy = dy + half_box_y end
         if (dz <-half_box_z) then dz = dz + half_box_z end
         var cutoff2 = subset[redo][part1].core_part_space.cutoff
         cutoff2 = cutoff2 * cutoff2
         var r2 = dx*dx + dy*dy + dz*dz
         regentlib.assert(r2 > 0.0, "Distance of 0 between particles")
         if(r2 <= cutoff2) then
           [nonsym_density_kernel(rexpr subset[redo][part1] end, rexpr parts[part2] end, rexpr r2 end)]
         end
       end
     end
   end
end

function reset_cutoff_update_space(part, space)
  local kernel = rquote
    var redo = int1d(0)
    part.cutoff_update_space.redo = redo
    part.cutoff_update_space.left = 0.0
    part.cutoff_update_space.right = hmax
    part.cutoff_update_space.h_0 = part.h
  end
  return kernel
end

function prepare_for_force(part, space)
  local kernel = rquote
    --Check if redo flag still active
    regentlib.assert(part.cutoff_update_space.redo == int1d(0), 
                     "A particle still needs redoing after redo stage done")
    --TODO Set timestep
    --Prepare for the force step
    var rho_inv = 1.0 / part.rho
    var h_inv = 1.0 / part.h
    var curl_v = sqrt( part.rot_v_x * part.rot_v_x +
                       part.rot_v_y * part.rot_v_y +
                       part.rot_v_z * part.rot_v_z )
    var div_v = part.div_v
    var abs_div_v = fabsd(div_v)
    var pressure = gas_pressure_from_internal_energy( part.rho, part.u )
    var soundspeed = gas_soundspeed_from_pressure( part.rho, pressure )
    var rho_dh = part.rho_dh
    --Ignore kernel changing affects if h is ~= hmax
    if(part.h > 0.9999 * hmax) then
      rho_dh = 0.0
    end 
  
    var grad_h_term = 1.0 / (1.0 + hydro_dimension_inv * part.h * rho_dh * rho_inv)
    var balsara = alpha * abs_div_v / (abs_div_v * curl_v + 0.0001 * soundspeed * h_inv)
    part.f = grad_h_term
    part.pressure = pressure
    part.soundspeed = soundspeed
    part.balsara = balsara
    regentlib.assert(1==0, "Updated h but not yet cutoff radius still")
    --Reset acceleration
    part.accel_x = 0.0
    part.accel_y = 0.0
    part.accel_z = 0.0
    part.u_dt = 0.0
    part.h_dt = 0.0
    part.v_sig = 2.0 * soundspeed
  end
  return kernel
end

function finish_density(part, space)
  local kernel = rquote
    var no_redo = int1d(1)
    var redo = int1d(0)
    var h_init = part.cutoff_update_space.h_0
    var h_old = part.h
    var h_old_dim = cube_val(h_old)
    var h_old_dim_minus_one = h_old * h_old
    var h_new : double
    var has_no_neighbours : bool
    var hydro_eta_3 = cube_val(hydro_eta)
    --Mark particles as done unless proven otherwise
    part.cutoff_update_space.redo = no_redo


    if(part.wcount == 0.0 ) then
      has_no_neighbours = true
      h_new = 2.0 * h_old
      regentlib.assert( 1 == 0, "Not yet implemented behaviour for neighbourless particles")
    else
      --Finish the density calculation
      var inv_h = 1.0 / part.h
      var inv_h_3 = cube_val(inv_h)
      var inv_h_4 = inv_h * inv_h_3
      --Add self contribution
      part.rho = part.rho + (part.core_part_space.mass * kernel_root)
      part.rho_dh = part.rho_dh - (kernel_dimension * part.core_part_space.mass * kernel_root)
      part.wcount = part.wcount + kernel_root
      part.wcount_dh = part.wcount_dh - (kernel_dimension * kernel_root)
      --Insert the missing h factor
      part.rho = part.rho * inv_h_3
      part.rho_dh = part.rho_dh * inv_h_4
      part.wcount = part.wcount * inv_h_3
      part.wcount_dh = part.wcount_dh - (kernel_dimension * kernel_root)
      --Complete the curl calculation
     var inv_rho = 1.0 / part.rho
     part.rot_v_x = part.rot_v_x * (inv_h_4 * inv_rho)
     part.rot_v_y = part.rot_v_y * (inv_h_4 * inv_rho)
     part.rot_v_z = part.rot_v_z * (inv_h_4 * inv_rho)
     part.div_v = part.div_v * (inv_h_4 * inv_rho)
 
     --Compute a Newton-Raphson step on h
     var n_sum = part.wcount * h_old_dim
     var n_target = hydro_eta_3
     var f = n_sum - n_target
     var f_prime = part.wcount_dh * h_old_dim + kernel_dimension * part.wcount * h_old_dim_minus_one
     --Improve the bounds
     if(n_sum < n_target) then
       part.cutoff_update_space.left = fmaxd(part.cutoff_update_space.left, h_old)
     elseif(n_sum > n_target) then
       part.cutoff_update_space.right = fmind(part.cutoff_update_space.right, h_old)
     end
     if( ( part.h >= hmax and f < 0.0 ) or (part.h <= hmin and f> 0.0) ) then
       regentlib.assert(1 == 0, "NYI: This should never happen as hmin is 0 and hmax is \"large\"")
      --FUTURE TODO: Calculate timestep and prepare for the force step and this particle is done.
     end
     --Finish the Newton-Raphson step on h
     h_new = h_old - f / (f_prime + 1e-128)
     --Limit the change to a factor of 2
     h_new = fmind(h_new, 2.0 * h_old)
     h_new = fmaxd(h_new, 0.5 * h_old)
     --Check for progression
     h_new = fmaxd(h_new, part.cutoff_update_space.left)
     h_new = fmind(h_new, part.cutoff_update_space.right)
    end

    --Now check whether the particle has an inappropriate smoothing length
    --If so, fix this and add the particle to the redo partition and reset the density parameters
    if( fabsd(h_new - h_old) > hydro_eps * h_old) then
      --Bisect if oscillating
      if( (h_new ==part.cutoff_update_space.left and h_old == part.cutoff_update_space.right) or (h_old == part.cutoff_update_space.left and h_new == part.cutoff_update_space.right)) then
        part.h = cbrtd( 0.5 * (cube_val(part.cutoff_update_space.left) + cube_val(part.cutoff_update_space.right) ) )
      else --Not oscillating
        part.h = h_new
      end --End of oscillation check
      --If the h is in an allowed range, then redo this part
      if(part.h < hmax and part.h > hmin) then
        --Make sure we know this needs redoing
        part.cutoff_update_space.redo = redo
        --Reset the properties we need for the density calculation
        part.rho = 0.0
        part.wcount = 0.0
        part.wcount_dh = 0.0
        part.rho_dh = 0.0
        part.div_v = 0.0
        part.rot_v_x = 0.0
        part.rot_v_y = 0.0
        part.rot_v_z = 00
        --TODO: Fix cutoff from H
        regentlib.assert(1 == 0, "Updated h but not yet cutoff radius")
      else
        --TODO: NYI as all particles should do previous check with current values
        regentlib.assert(1==0, "NYI: Particles outside acceptable range")
      end
    end

  end
  --At the end of this kernel, particles redo values are either:
  --1) redo, with 0d density values ready to be redone
  --2) no_redo, ready to have the values for the force loop computed.
  return kernel
end


local reset_cutoff_update_task_runner = run_per_particle_task( reset_cutoff_update_space )
local prepare_for_force_runner = run_per_particle_task(  prepare_for_force )

local sort_redo_task = generate_per_part_task( finish_density )

--The SPH task between density and force calculations is a special task, which for now is 
--non-generalised.
task update_cutoffs_launcher(particles : region(ispace(int1d), part),
                             cell_space : partition(disjoint, particles, ispace(int3d)),
                             space : region(ispace(int1d), space_config)) where
                             reads(particles, space), writes(particles) do
    var no_redo = int1d(1)
    var redo = int1d(0)

    --First we reset all the particle's cutoffs
    reset_cutoff_update_task_runner( particles, cell_space, space)
    --Loop over particles with redo set until all particles are redone
    for attempts = 1, 5 do
      var redo_partition = partition(particles.cutoff_update_space.redo, ispace(int1d, 2))
      --We have two partitions, redo partition and the cell partition, we want the cross product of this.
      var cell_redo = cross_product(cell_space, redo_partition)
      for cell in cell_space.colors do
        sort_redo_task(cell_redo[cell][redo], space)
      end
      __delete(redo_partition)
    end

    --At the end of the loop everyone should be redone succesfully
    prepare_for_force_runner(particles, cell_space, space)
end


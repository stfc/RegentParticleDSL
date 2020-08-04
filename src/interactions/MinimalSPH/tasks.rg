import "regent"

require("defaults")
--require("src/neighbour_search/cell_pair/import_cell_pair")
require("src/neighbour_search/cell_pair/neighbour_search")
require("src/neighbour_search/cell_pair/cell")
require("src/interactions/MinimalSPH/interactions")
require("src/interactions/MinimalSPH/timestep")

local density_symmetric_task = generate_asymmetric_pairwise_task(density_kernel)
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

--The SPH task between density and force calculations is a special task, which for now is 
--non-generalised.
task update_cutoffs(parts1 : region(ispace(int1d),part), particles: region(ispace(int1d), part), cell_space : partition(disjoint, particles , ispace(int3d)), space : region(ispace(int1d), space_config), this_cell : int3d ) where reads(parts1, particles, space), writes(parts1) do
--Define some constants/import some math functions
var no_redo = int1d(1)
var redo = int1d(0)
var hydro_eta_3 = cube_val(hydro_eta)

for particle in parts1 do
  --Reset values
  parts1[particle].cutoff_update_space.redo = redo
  parts1[particle].cutoff_update_space.left = 0.0
  parts1[particle].cutoff_update_space.right = hmax
  parts1[particle].cutoff_update_space.h_0 = parts1[particle].h
end


--TODO: Loop more than once
for attempts = 1, 2 do
  var redo_partition = partition(parts1.cutoff_update_space.redo, ispace(int1d, 2))
  --Loop through the particles and correct h
  for particle in redo_partition[redo] do
    var h_init = parts1[particle].cutoff_update_space.h_0
    var h_old = parts1[particle].h
    var h_old_dim = cube_val(h_old)
    var h_old_dim_minus_one = h_old * h_old
    var h_new : double
    var has_no_neighbours : bool
    parts1[particle].cutoff_update_space.redo = no_redo

    if( parts1[particle].wcount == 0.0 ) then
      has_no_neighbours = true
      --Double the cutoff and retry
      h_new = 2.0 * h_old
      regentlib.assert(1 == 0, "Not yet implementeds stuff happens here")
    else
      --TODO: Finish the density calculation
      var inv_h = 1.0 / parts1[particle].h
      var inv_h_3 = cube_val(inv_h)
      var inv_h_4 = inv_h * inv_h_3
   
      --Add self contribution
      parts1[particle].rho = parts1[particle].rho + ( parts1[particle].core_part_space.mass * kernel_root)
      parts1[particle].rho_dh = parts1[particle].rho_dh - ( kernel_dimension * parts1[particle].core_part_space.mass * kernel_root )
      parts1[particle].wcount = parts1[particle].wcount + kernel_root
      parts1[particle].wcount_dh = parts1[particle].wcount_dh - ( kernel_dimension * kernel_root )
      --Insert the missing h factors
      parts1[particle].rho = parts1[particle].rho * inv_h_3
      parts1[particle].rho_dh = parts1[particle].rho_dh * inv_h_4
      parts1[particle].wcount = parts1[particle].wcount * inv_h_3
      parts1[particle].wcount_dh = parts1[particle].wcount_dh * inv_h_4 
      --Complete the curl calculation
      var inv_rho = 1.0 / parts1[particle].rho
      parts1[particle].rot_v_x = parts1[particle].rot_v_x * ( inv_h_4 * inv_rho)
      parts1[particle].rot_v_y = parts1[particle].rot_v_y * ( inv_h_4 * inv_rho)
      parts1[particle].rot_v_z = parts1[particle].rot_v_z * ( inv_h_4 * inv_rho)
      parts1[particle].div_v = parts1[particle].div_v * ( inv_h_4 * inv_rho)

      --TODO: Compute a newton-raphson step on h
      var n_sum = parts1[particle].wcount * h_old_dim
      var n_target = hydro_eta_3
      var f = n_sum - n_target
      var f_prime = parts1[particle].wcount_dh * h_old_dim + kernel_dimension * parts1[particle].wcount * h_old_dim_minus_one
        --Improve the bounds
        if(n_sum < n_target) then
          parts1[particle].cutoff_update_space.left = fmaxd(parts1[particle].cutoff_update_space.left, h_old)
        elseif(n_sum > n_target) then
          parts1[particle].cutoff_update_space.right = fmind(parts1[particle].cutoff_update_space.right, h_old)
        end
      --TODO: Finish with a particle if above hmax or below h_min sometimes
      if( (parts1[particle].h >= hmax and f < 0.0 ) or (parts1[particle].h <= hmin and f > 0.0) ) then
        regentlib.assert(1 == 0, "NYI: This should never happen as hmin is 0 and hmax is \"large\"")
        --FUTURE TODO: Calculate timestep and prepare for the force step and this particle is done.
      end

      --TODO: Finish the newton-raphson step on h
          h_new = h_old - f / (f_prime + 1e-128);
         --Limit the change to a factor of 2
         h_new = fmind(h_new, 2.0 * h_old)
         h_new = fmaxd(h_new, 0.5 * h_old)
         --Check for progression
         h_new = fmaxd(h_new, parts1[particle].cutoff_update_space.left)
         h_new = fmind(h_new, parts1[particle].cutoff_update_space.right)

    end
      --TODO: Check whether the particle has an inappropriate smoothing length
      --If so fix this and add the particle to the redo partition. Reset the density parameters
      if( fabsd(h_new - h_old) > hydro_eps * h_old) then
        --Bisect if we're oscillating
        if( (h_new == parts1[particle].cutoff_update_space.left and h_old == parts1[particle].cutoff_update_space.right) or
             (h_old == parts1[particle].cutoff_update_space.left and h_new == parts1[particle].cutoff_update_space.right) ) then
           parts1[particle].h = cbrtd( 0.5 * (cube_val(parts1[particle].cutoff_update_space.left) + cube_val(parts1[particle].cutoff_update_space.right) ) )
        else
           parts1[particle].h = h_new
        end
        --If in acceptable range then redo this part
        if( parts1[particle].h <hmax and parts1[particle].h > hmin) then
          parts1[particle].cutoff_update_space.redo = redo
          --Reset the particle properties
          parts1[particle].rho = 0.0
          parts1[particle].wcount = 0.0
          parts1[particle].wcount_dh = 0.0
          parts1[particle].rho_dh = 0.0
          parts1[particle].div_v = 0.0
          parts1[particle].rot_v_x = 0.0
          parts1[particle].rot_v_y = 0.0
          parts1[particle].rot_v_z = 0.0
          --TODO FIX CUTOFF FROM H
          regentlib.assert(1 == 0, "Updated h but not yet cutoff radius")
          else
            --TODO: NYI as all particles should do previous check for now
            regentlib.assert(1==0, "NYI: Particles outside acceptable range")
        end
      end

  end

  --Recreate the redo partition
  __delete(redo_partition)
  var redo_partition2 = partition(parts1.cutoff_update_space.redo, ispace(int1d, 2))  
      --TODO: For particles that have a converged smoothing length, get ready for force loop.
  for particle in redo_partition2[no_redo] do
    --TODO Set timestep
    --Prepare force
    var rho_inv = 1.0 / redo_partition2[no_redo][particle].rho
    var h_inv = 1.0 / redo_partition2[no_redo][particle].h
    var curl_v = sqrt(redo_partition2[no_redo][particle].rot_v_x * redo_partition2[no_redo][particle].rot_v_x + 
                      redo_partition2[no_redo][particle].rot_v_y * redo_partition2[no_redo][particle].rot_v_y + 
                      redo_partition2[no_redo][particle].rot_v_z * redo_partition2[no_redo][particle].rot_v_z)
    var div_v = redo_partition2[no_redo][particle].div_v
    var abs_div_v = fabsd(div_v)
    var pressure = gas_pressure_from_internal_energy( redo_partition2[no_redo][particle].rho, redo_partition2[no_redo][particle].u )
    var soundspeed = gas_soundspeed_from_pressure(  redo_partition2[no_redo][particle].rho, pressure )
    var rho_dh = redo_partition2[no_redo][particle].rho_dh
    --Ignore kernel changing effects if h approx equal hmax
    if(redo_partition2[no_redo][particle].h > 0.9999 * hmax) then
       rho_dh = 0.0
    end
    
    var grad_h_term = 1.0 / (1.0 + hydro_dimension_inv * redo_partition2[no_redo][particle].h * rho_dh * rho_inv)
    var balsara = alpha * abs_div_v / ( abs_div_v * curl_v + 0.0001 * soundspeed * h_inv)
    redo_partition2[no_redo][particle].f = grad_h_term
    redo_partition2[no_redo][particle].pressure = pressure
    redo_partition2[no_redo][particle].soundspeed = soundspeed
    redo_partition2[no_redo][particle].balsara = balsara 

          regentlib.assert(1 == 0, "Updated h but not yet cutoff radius")
    --TODO: Reset acceleration
    redo_partition2[no_redo][particle].accel_x = 0.0
    redo_partition2[no_redo][particle].accel_y = 0.0
    redo_partition2[no_redo][particle].accel_z = 0.0
    redo_partition2[no_redo][particle].u_dt = 0.0
    redo_partition2[no_redo][particle].h_dt = 0.0
    redo_partition2[no_redo][particle].v_sig = 2.0 * soundspeed
  end
  --TODO: Loop over the redo particles and rerun the density loops on this subset with all other cells.
  --FUTURE TODO: Can use a generated task to do this for an appropriate neighbour search algorithm 
--  run_subset_task(redo_partition2[redo], particles, cell_space, space )
  --Wait on the subset task -- TODO: May not need this
--  self_redo_density( parts1, redo_partition2, space )
  for cells in cell_space.colors do
    if cells ~= this_cell then
      pair_redo_density( parts1 , redo_partition2, cell_space[cells], space )
    end
  end
  c.legion_runtime_issue_execution_fence(__runtime(), __context())
  __delete(redo_partition2)
end

end

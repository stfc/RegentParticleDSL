import "regent"

require("defaults")
require("src/neighbour_search/cell_pair/neighbour_search")
require("src/neighbour_search/cell_pair/cell")
require("src/interactions/MinimalSPH/interactions")
require("src/interactions/MinimalSPH/timestep")

local density_task = create_asymmetric_pairwise_runner(nonsym_density_kernel)
local c = regentlib.c
local stdlib = terralib.includec("stdlib.h")

local ceil = regentlib.ceil(float)

--Some big number
local hmax = 12345678.0
local hmin = 0.0

local fabsd = regentlib.fabs(float)
local cbrtd = regentlib.cbrt(float)
local sqrt = regentlib.sqrt(float)
--3D case
terra cube_val(value : float)
  var x2 = value * value
  return x2 * value
end

terra fmind(val1 : float, val2 : float)
  return regentlib.fmin(val1, val2)
end

terra fmaxd(val1 : float, val2 : float)
  return regentlib.fmax(val1, val2)
end

terra allocate_bool_array( count : int ) 
  return stdlib.malloc(count * sizeof(bool) )
end

local terra create_temp_array(size : int)
  return stdlib.malloc(size * sizeof(double))
end

local terra free_array(array : &double)
  return stdlib.free(array)
end
                       
task compute_timesteps(particles: region(ispace(int1d), part)) : double where
  reads(particles.h, particles.v_sig) do
  var min_timestep = 1000000000.0
  for part in particles.ispace do
   var dt_cfl = 2.0 * kernel_gamma * CFL_condition * particles[part].h / ( particles[part].v_sig)
    min_timestep = regentlib.fmin(min_timestep, dt_cfl)
  end
  return min_timestep
end

task compute_timestep_launcher( particles: region(ispace(int1d), part), cell_space : partition(disjoint, particles , ispace(int3d)), config : region(ispace(int1d), config_type) ) : double
    where reads(particles.h, particles.v_sig) do
   
   var dx = cell_space.colors.bounds.hi.x - cell_space.colors.bounds.lo.x + 1
   var dy = cell_space.colors.bounds.hi.y - cell_space.colors.bounds.lo.y + 1
   var dz = cell_space.colors.bounds.hi.z - cell_space.colors.bounds.lo.z + 1

   var count = dx*dy*dz
   var array : &double = [&double](create_temp_array(count))

    count = 0
    --For each cell, call the task!
    for cell1 in cell_space.colors do
       array[count] = compute_timesteps(cell_space[cell1])
       count = count + 1
    end

     count = 0
    var min_timestep = 1000000000.0
    --For each cell, call the task!
    for cell1 in cell_space.colors do
       min_timestep = regentlib.fmin(min_timestep, array[count])
       count = count + 1
    end

    free_array(array)
    return min_timestep
end

task pair_redo_density( parts_self : region(ispace(int1d), part), 
                        parts_far : region(ispace(int1d), part),
                        config : region(ispace(int1d), config_type) )
   where reads(parts_self, parts_far, config), writes( parts_self ) do
   var no_redo = int1d(1)
   var redo = int1d(0)
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var box_z = config[0].space.dim_z
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   var half_box_z = 0.5 * box_z
   for part1 in parts_self.ispace do
     if(parts_self[part1].cutoff_update_space.redo == redo) then
       for part2 in parts_far.ispace do
         --Compute particle distance
           var dx = parts_self[part1].core_part_space.pos_x - parts_far[part2].core_part_space.pos_x
           var dy = parts_self[part1].core_part_space.pos_y - parts_far[part2].core_part_space.pos_y
           var dz = parts_self[part1].core_part_space.pos_z - parts_far[part2].core_part_space.pos_z
           if (dx > half_box_x) then dx = dx - box_x end
           if (dy > half_box_y) then dy = dy - box_y end
           if (dz > half_box_z) then dz = dz - box_z end
           if (dx <-half_box_x) then dx = dx + box_x end
           if (dy <-half_box_y) then dy = dy + box_y end
           if (dz <-half_box_z) then dz = dz + box_z end
           var cutoff2 = parts_self[part1].core_part_space.cutoff
           cutoff2 = cutoff2 * cutoff2
           var r2 = dx*dx + dy*dy + dz*dz
           if(r2 <= cutoff2) then
             [nonsym_density_kernel(rexpr parts_self[part1] end, rexpr parts_far[part2] end, rexpr r2 end)]
           end
       end
     end
   end
end

task self_redo_density( parts : region(ispace(int1d), part), config : region(ispace(int1d), config_type) )
   where reads(parts, config), writes(parts) do
   var no_redo = int1d(1)
   var redo = int1d(0)
   var box_x = config[0].space.dim_x
   var box_y = config[0].space.dim_y
   var box_z = config[0].space.dim_z
   var half_box_x = 0.5 * box_x
   var half_box_y = 0.5 * box_y
   var half_box_z = 0.5 * box_z
   for part1 in parts.ispace do
     if(parts[part1].cutoff_update_space.redo == redo) then
       for part2 in parts.ispace do
         --Compute particle distance
         if(int1d(part1) ~= int1d(part2)) then
           var dx = parts[part1].core_part_space.pos_x - parts[part2].core_part_space.pos_x
           var dy = parts[part1].core_part_space.pos_y - parts[part2].core_part_space.pos_y
           var dz = parts[part1].core_part_space.pos_z - parts[part2].core_part_space.pos_z
           if (dx > half_box_x) then dx = dx - box_x end
           if (dy > half_box_y) then dy = dy - box_y end
           if (dz > half_box_z) then dz = dz - box_z end
           if (dx <-half_box_x) then dx = dx + box_x end
           if (dy <-half_box_y) then dy = dy + box_y end
           if (dz <-half_box_z) then dz = dz + box_z end
           var cutoff2 = parts[part1].core_part_space.cutoff
           cutoff2 = cutoff2 * cutoff2
           var r2 = dx*dx + dy*dy + dz*dz
           regentlib.assert(r2 > 0.0, "Distance of 0 between particles")
           if(r2 <= cutoff2) then
             [nonsym_density_kernel(rexpr parts[part1] end, rexpr parts[part2] end, rexpr r2 end)]
           end
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
    var string : rawstring
    string = [rawstring](stdlib.malloc( 512 )) 
    format.snprintln(string, 512, "A particle still needs redoing after redo stage done id = {}", part.core_part_space.id)
    regentlib.assert(part.cutoff_update_space.redo == int1d(1), 
                     string)
    stdlib.free(string)
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
    part.core_part_space.cutoff = kernel_gamma * part.h
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

function finish_density(part, space, return_bool)
  local kernel = rquote
  var no_redo = int1d(1)
  var redo = int1d(0)
  if(part.cutoff_update_space.redo == redo) then 
    var h_init : float = part.cutoff_update_space.h_0
    var h_old : float = part.h
    var h_old_dim : float = cube_val(h_old)
    var h_old_dim_minus_one : float = h_old * h_old
    var h_new : float
    var has_no_neighbours : bool
    var hydro_eta_3 = cube_val(hydro_eta)
    --Mark particles as done unless proven otherwise
    part.cutoff_update_space.redo = no_redo

    if(part.wcount == 0.0 ) then
      has_no_neighbours = true
      h_new = 2.0 * h_old
      var string : rawstring
      string = [rawstring](stdlib.malloc( 512 )) 
      format.snprintln(string, 512, "neighbourless particles wcount = {} cutoff  = {} h_init = {} id = {}", part.wcount, h_old, part.cutoff_update_space.h_0, part.core_part_space.id)
      regentlib.assert( 1 == 0, string )--"Not yet implemented behaviour for neighbourless particles")
      stdlib.free(string)
    else
      --Finish the density calculation
      var inv_h : float = 1.0 / part.h
      var inv_h_3 : float = cube_val(inv_h)
      var inv_h_4 : float = inv_h * inv_h_3
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
        part.core_part_space.cutoff = kernel_gamma * part.h
      else --Not oscillating
        part.h = h_new
        part.core_part_space.cutoff = kernel_gamma * part.h
      end --End of oscillation check
      --If the h is in an allowed range, then redo this part
      if(part.h < hmax and part.h > hmin) then
        --Make sure we know this needs redoing
        part.cutoff_update_space.redo = redo
        return_bool = true
        --Reset the properties we need for the density calculation
        part.rho = 0.0
        part.wcount = 0.0
        part.wcount_dh = 0.0
        part.rho_dh = 0.0
        part.div_v = 0.0
        part.rot_v_x = 0.0
        part.rot_v_y = 0.0
        part.rot_v_z = 00
        part.core_part_space.cutoff = kernel_gamma * part.h
      else
        --TODO: NYI as all particles should do previous check with current values
        regentlib.assert(1==0, "NYI: Particles outside acceptable range")
      end
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

local sort_redo_task = generate_per_part_task_bool_return( finish_density )

--The SPH task between density and force calculations is a special task, which for now is 
--non-generalised.
task update_cutoffs_launcher(particles : region(ispace(int1d), part),
                             cell_space : partition(disjoint, particles, ispace(int3d)),
                             config : region(ispace(int1d), config_type)) where
                             reads(particles, config), writes(particles) do
    var no_redo = int1d(1)
    var redo = int1d(0)
    --Compute cell radii
    var x_count = config[0].neighbour_config.x_cells
    var y_count = config[0].neighbour_config.y_cells
    var z_count = config[0].neighbour_config.z_cells
    var cutoff = config[0].neighbour_config.max_cutoff
    --Hacky fix to deal with increasing cutoffs
    var x_radii : int = ceil( cutoff / config[0].neighbour_config.cell_dim_x ) +1
    var y_radii : int = ceil( cutoff / config[0].neighbour_config.cell_dim_y ) +1
    var z_radii : int = ceil( cutoff / config[0].neighbour_config.cell_dim_z ) +1
    var nr_cells = x_count * y_count * z_count
    var bool_array : &bool = [&bool](allocate_bool_array(nr_cells) )
    for i=0,nr_cells do
      bool_array[i] = true
    end
    for cell in cell_space.colors do
      var point = cell.x + cell.y * x_count + cell.z * x_count * y_count
    end
    --First we reset all the particle's cutoffs
    reset_cutoff_update_task_runner( particles, cell_space, config)
    --Loop over particles with redo set until all particles are redone
    for attempts = 1, 10 do
      var launches = 0
      for cell in cell_space.colors do
        var point = cell.x + cell.y * x_count + cell.z * x_count * y_count
        if(bool_array[point]) then
          bool_array[point] = sort_redo_task(cell_space[cell], config)
          launches = launches + 1
        end
      end
c.legion_runtime_issue_execution_fence(__runtime(), __context())
      for cell1 in cell_space.colors do
        var point = cell1.x + cell1.y * x_count + cell1.z * x_count * y_count
        if(bool_array[point]) then
          self_redo_density(cell_space[cell1], config)
        end
        for x = 0, x_radii+1 do
          for y = 0, y_radii+1 do
            for z = 0, z_radii + 1 do
              if(not (x == 0 and y == 0 and z == 0) ) then
                var cell2 : int3d = int3d({ (cell1.x + x)%x_count, (cell1.y +y)%y_count, (cell1.z + z)%z_count })
                var point2 = cell2.x + cell2.y * x_count + cell2.z * x_count * y_count
                --Weird if statement to handle max_cutoff >= half the boxsize
                if( (cell1.x > cell2.x or (cell1.x == cell2.x and cell1.y > cell2.y) or( cell1.x == cell2.x and cell1.y == cell2.y and cell1.z > cell2.z)) and
                   (cell1.x - x_radii <= cell2.x and cell1.y - y_radii <= cell2.y and cell1.z - y_radii <= cell2.z) ) then
           
                else
                  --Asymmetric, launch tasks both ways
                  if(bool_array[point]) then
                    pair_redo_density( cell_space[cell1], cell_space[cell2], config)
                  end
                  if(bool_array[point2]) then
                    pair_redo_density( cell_space[cell2], cell_space[cell1], config)
                  end
                end
              end
            end
          end
        end
      end
    end
      for cell in cell_space.colors do
        var point = cell.x + cell.y * x_count + cell.z * x_count * y_count
        if(bool_array[point]) then
          bool_array[point] = sort_redo_task(cell_space[cell], config)
        end
      end
c.legion_runtime_issue_execution_fence(__runtime(), __context())
    stdlib.free(bool_array)
    --At the end of the loop everyone should be redone succesfully
    prepare_for_force_runner(particles, cell_space, config)
end


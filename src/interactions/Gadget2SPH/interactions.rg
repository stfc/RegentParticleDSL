import "regent"


require("defaults")
require("src/neighbour_search/cell_pair/neighbour_search")
local sqrt = regentlib.sqrt(double)



--TODO: Optimise kernel evaluation (Could use C code?)
--Think this should be cubic spline, 3 Dimension form only.
terra eval_kernel_func( ui : double, wi : double[2] )
  var kernel_gamma_inv = 1.0 / 1.825742
  var q3 = 0.0
  var q2 = 0.0
  var q1 = 0.0
  var q0 = 0.0
  var x = ui * kernel_gamma_inv

  if ui <= 0.5 then
    q3 = 0.75
    q2 = 1.5
    q1 = 0
    q0 = 1.0
  else
    q3 = -0.25
    q2 = 1.5
    q1 = -3.0
    q0 = 2
  end

  wi[0] = q3 * x + q2
  wi[1] = q3
  --
  wi[1] = wi[1] * x + wi[0]
  wi[0] = wi[0] * x + q1
  --
  wi[1] = wi[1] * x + wi[0]
  wi[0] = wi[0] * x + q0
end


function iact_density(part1, part2, r2)

local hydro_dimension = 3.0
local kernel = rquote
  var mi = part1.mass
  var mj = part2.mass

  var ir = 1.0 / sqrt(r2)
  var r = r2 * ir

  var hi_inv = 1.0 / part1.core_part_space.cutoff
  var ui = r * hi_inv
  var wi_array : double[2]
  eval_kernel_func(ui, wi_array)
  var wi : double = wi_array[0]
  var wi_dx : double = wi_array[1]

  part1.rho = (part1.rho + mj * wi)
  part1.rho_dh = part1.rho_dh - (mj * hydro_dimension * wi + ui * wi_dx)
  part1.wcount = part1.wcount + wi
  part1.wcount_dh = part1.wcount_dh - (hydro_dimension * wi + ui * wi_dx)

  var hj_inv = 1.0 / part2.core_part_space.cutoff
  var uj = r * hj_inv
  var wj_array : double[2]
  eval_kernel_func(uj, wj_array)
  var wj : double = wj_array[0]
  var wj_dx : double = wj_array[1]
  part2.rho = part2.rho + mi * wj
  part2.rho_dh = part2.rho_dh - (mi * hydro_dimension * wj + uj * wj_dx)
  part2.wcount = part2.wcount + wj
  part2.wcount_dh = part2.wcount_dh - (hydro_dimension * wj + uj * wj_dx)

  var faci = mj * wi_dx * ir
  var facj = mi * wj_dx * ir

  --Compute dvdr
  var dv_x = part1.vel_x - part2.vel_x
  var dv_y = part1.vel_y - part2.vel_y
  var dv_z = part1.vel_z - part2.vel_z
  var dx_x = part1.core_part_space.pos_x - part2.core_part_space.pos_x
  var dx_y = part1.core_part_space.pos_y - part2.core_part_space.pos_y
  var dx_z = part1.core_part_space.pos_z - part2.core_part_space.pos_z
  var dvdr = dv_x * dx_x + dv_y * dx_y + dv_z * dx_z

  part1.div_v = part1.div_v - (faci * dvdr)
  part2.div_v = part2.div_v - (facj * dvdr)

  --Compute dv cross r
  var curlvr_x = dv_y * dx_z - dv_z * dx_y
  var curlvr_y = dv_z * dx_x - dv_x * dx_z
  var curlvr_z = dv_x * dx_y - dv_y * dx_x

  part1.rot_v_x = part1.rot_v_x + (faci * curlvr_x)
  part1.rot_v_y = part1.rot_v_y + (faci * curlvr_y)
  part1.rot_v_z = part1.rot_v_z + (faci * curlvr_z)
  
 
  part2.rot_v_x = part2.rot_v_x + (facj * curlvr_x)
  part2.rot_v_y = part2.rot_v_y + (facj * curlvr_y)
  part2.rot_v_z = part2.rot_v_z + (facj * curlvr_z)
end
return kernel
end

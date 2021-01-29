import "regent"

local sqrtf = regentlib.sqrt(float)


--TODO Implement 2D version

terra pow_dimension_plus_one(x : float) : float
  var x2 = x*x
  return x2*x2
end


local M_PI = 3.141592653589793238462643383279502884
local kernel_gamma = global(float, 1.825742)
local kernel_constant = global(float, 16.0 / M_PI)
local kernel_gamma_inv_dim = global(float, 1.0 / (kernel_gamma:get() * kernel_gamma:get() * kernel_gamma:get()) )
local kernel_gamma_inv_dim_plus_one = global(float, 1.0 / (kernel_gamma:get() * kernel_gamma:get() * kernel_gamma:get() * kernel_gamma:get()) )

terra kernel_deval( ui : float, wi : &float, wi_dx : &float )
  var kernel_gamma_inv = 1.0 / kernel_gamma
  var q3 = 0.0
  var q2 = 0.0
  var q1 = 0.0
  var q0 = 0.0
  var x = ui * kernel_gamma_inv
  var w :float = 0.0
  var w_dx : float = 0.0
  if ui <= 0.5 then
    q3 = 3.
    q2 = -3.0
    q1 = 0
    q0 = 0.5
  else
    q3 = -1.0
    q2 = 3.0
    q1 = -3.0
    q0 = 1.0
  end

  w = q3 * x + q2
  w_dx = q3

  w_dx = w_dx * x + w
  w = w * x + q1

  w_dx = w_dx * x + w
  w = w * x + q0

  if (w < 0.0) then w = 0.0 end
  if (w_dx > 0.0) then w_dx = 0.0 end
  w = w * kernel_constant * kernel_gamma_inv_dim
  w_dx = w_dx * kernel_constant * kernel_gamma_inv_dim_plus_one

  @wi = w
  @wi_dx = w_dx
end

--Default: Laminar SPS Kernel
function force_kernel(part1, part2, r2)
local kernel = rquote
  var r : float = sqrtf(r2)
  var mi : float = part1.core_part_space.mass
  var mj : float = part2.core_part_space.mass
  var rhoi : float = part1.rho
  var rhoj : float = part2.rho
  var rhoi_inv : float = 1.0 / rhoi
  var rhoj_inv : float = 1.0 / rhoj
  var r_inv : float = 1.0 / r

  var dx0 : float = part1.core_part_space.pos_x - part2.core_part_space.pos_x
  var dx1 : float = part1.core_part_space.pos_y - part2.core_part_space.pos_y
  var dx2 : float = part1.core_part_space.pos_z - part2.core_part_space.pos_z
  var dv0 : float = part1.core_part_space.vel_x - part2.core_part_space.vel_x
  var dv1 : float = part1.core_part_space.vel_y - part2.core_part_space.vel_y
  var dv2 : float = part1.core_part_space.vel_z - part2.core_part_space.vel_z

  var inv_hidim_pow_plus_one : float = 1.0 / pow_dimension_plus_one(part1.h)
  var inv_hjdim_pow_plus_one : float = 1.0 / pow_dimension_plus_one(part2.h)

  var hi_inv : float = 1.0 / part1.h
  var xi : float= r * hi_inv
  var wi : float = 0.0
  var wi_dx : float = 0.0
  kernel_deval(xi, &wi, &wi_dx)
  wi_dx = wi_dx  * inv_hidim_pow_plus_one

  var hj_inv : float = 1.0 / part2.h
  var xj : float= r * hj_inv
  var wj : float = 0.0
  var wj_dx : float = 0.0
  kernel_deval(xj, &wj, &wj_dx)
  wj_dx = wj_dx  * inv_hjdim_pow_plus_one

  var r2_eta2 : float = r2 + 0.01*0.01
  var inv_r2eta2 : float = 1.0 / r2_eta2
  var visc_i : float = 4.0 * part1.viscosity
  var visc_j : float = 4.0 * part2.viscosity

  var temp_i : float = visc_i / ( (r2_eta2) * (rhoi + rhoj))
  var temp_j : float = visc_j / ( (r2_eta2) * (rhoi + rhoj))
  var multiplier_i : float = dx0*dx0*wi_dx + dx1*dx1*wi_dx + dx2*dx2*wi_dx
  var multiplier_j : float = dx0*dx0*wj_dx + dx1*dx1*wj_dx + dx2*dx2*wj_dx

  part1.a_hydro_x += (mj*temp_i*multiplier_i*dv0)
  part1.a_hydro_y += (mj*temp_i*multiplier_i*dv1)
  part1.a_hydro_z += (mj*temp_i*multiplier_i*dv2)

  part2.a_hydro_x -= (mi*temp_j*multiplier_j*dv0)
  part2.a_hydro_y -= (mi*temp_j*multiplier_j*dv1)
  part2.a_hydro_z -= (mi*temp_j*multiplier_j*dv2)

  --Compute turbulence
  var tau_xx : float = part1.tau_xx + part2.tau_xx
  var tau_xy : float = part1.tau_xy + part2.tau_xy
  var tau_xz : float = part1.tau_xz + part2.tau_xz
  var tau_yy : float = part1.tau_yy + part2.tau_yy
  var tau_yz : float = part1.tau_yz + part2.tau_yz
  var tau_zz : float = part1.tau_zz + part2.tau_zz
  var mi_mj : float = mi * mj
  var hydro0 : float =  mi_mj * (tau_xx*dx0 + tau_xy * dx1 + tau_xz * dx2);
  var hydro1 : float =  mi_mj * (tau_xy*dx0 + tau_yy * dx1 + tau_yz * dx2);
  var hydro2 : float =  mi_mj * (tau_xz*dx0 + tau_yz * dx1 + tau_zz * dx2);
  
  part1.a_hydro_x += hydro0
  part1.a_hydro_y += hydro1
  part1.a_hydro_z += hydro2

  part2.a_hydro_x -= hydro0
  part2.a_hydro_y -= hydro1
  part2.a_hydro_z -= hydro2

  --Compute velocity gradiants
  var mj_over_rhoj : float = mj * rhoj_inv
  var mi_over_rhoi : float = mi * rhoi_inv
  var dv_over_rho_x_i : float = dv0 * mj_over_rhoj
  var dv_over_rho_y_i : float = dv1 * mj_over_rhoj
  var dv_over_rho_z_i : float = dv2 * mj_over_rhoj
  var dv_over_rho_x_j : float = dv0 * mi_over_rhoi
  var dv_over_rho_y_j : float = dv1 * mi_over_rhoi
  var dv_over_rho_z_j : float = dv2 * mi_over_rhoi

  part1.grad_v_xx += (dv_over_rho_x_i * dx0)
  part1.grad_v_xy += (dv_over_rho_x_i * dx1)
  part1.grad_v_xz += (dv_over_rho_x_i * dx2)
  part1.grad_v_xy += (dv_over_rho_y_i * dx0)
  part1.grad_v_yy += (dv_over_rho_y_i * dx1)
  part1.grad_v_yz += (dv_over_rho_y_i * dx2)
  part1.grad_v_xz += (dv_over_rho_z_i * dx0)
  part1.grad_v_yz += (dv_over_rho_z_i * dx1)
  part1.grad_v_zz += (dv_over_rho_z_i * dx2)

  part2.grad_v_xx += (dv_over_rho_x_j * dx0)
  part2.grad_v_xy += (dv_over_rho_x_j * dx1)
  part2.grad_v_xz += (dv_over_rho_x_j * dx2)
  part2.grad_v_xy += (dv_over_rho_y_j * dx0)
  part2.grad_v_yy += (dv_over_rho_y_j * dx1)
  part2.grad_v_yz += (dv_over_rho_y_j * dx2)
  part2.grad_v_xz += (dv_over_rho_z_j * dx0)
  part2.grad_v_yz += (dv_over_rho_z_j * dx1)
  part2.grad_v_zz += (dv_over_rho_z_j * dx2)

  var acc : float = (part1.pressure + part2.pressure) * ((rhoi_inv * rhoj_inv) * r_inv)

  var dens : float = dv0*wi_dx*dx0 + dv1*wi_dx*dx1 + dv2*wi_dx*dx2

--  if(part1.core_part_space.id == int1d(1116)) then
--    format.println("hydro_x interaction {} {} {} {}", mj, acc, wi_dx, dx1)
--  end
  part1.drho_dt += (mj * dens)
  part1.a_hydro_x -= (mj * acc * wi_dx * dx0)
  part1.a_hydro_y -= (mj * acc * wi_dx * dx1)
  part1.a_hydro_z -= (mj * acc * wi_dx * dx2)

--  if(part2.core_part_space.id == int1d(1116)) then
--    format.println("hydro_y interaction {} {} {} {}", mi, acc, wj_dx, dx1)
--  end
  --if(part2.core_part_space.id == int1d(0)) then
  --  format.println("drho_dt interaction {} {}", mi, dens)
  --  var x = 0
  --  if  part1.neighbour_part_space._valid then
  --     x = 1
  --  end
  --  format.println("other part is real? {}", x)
  --end
  part2.drho_dt += (mi*dens)
  part2.a_hydro_x += (mi * acc * wj_dx * dx0)
  part2.a_hydro_y += (mi * acc * wj_dx * dx1)
  part2.a_hydro_z += (mi * acc * wj_dx * dx2)

  --Store viscous effect towards CFL condition
  var dvdr : float = dv0*dx0 + dv1*dx1 + dv2*dx2
  var dvdr_rr2 : float = dvdr * inv_r2eta2
  --FIXME : Sort out MAX function
  part1.max_visc max= max(part1.max_visc, dvdr_rr2)
  part2.max_visc max= max(part2.max_visc, dvdr_rr2)

  part1.interactions += 1
  part2.interactions += 1
  --FIXME: IMPLEMENT THESE
  [boundary_fluid_interaction( part1, part2, r, r2, dx0, dx1, dx2)];
  [compute_density_diffusive_term(part1, part2, r2, r, wi_dx, wj_dx, dx)];
  [compute_shifting_term(part1, part2, r2, r, wi_dx, wj_dx, dx)];
end
return kernel
end

import "regent"

require("defaults")

local sqrtf = regentlib.sqrt(float)


--TODO Implement 2D version

terra pow_dimension_plus_one(x : float) : float
  var x2 = x*x
  return x2*x2
end


--Cubic spline kernel
terra kernel_deval(u : float, h : float, w : &float, w_dx : &float)
  
  var u2 = u*u
  var u3 = u2*u
  var a_D = (10.0/7.0) * 3.141592654 * h*h
  
  if w >= 0 and w <= 1 then
    &w = a_D* 1.0 - 1.5*(u2) + 0.75*(u3)
    &w_dx = a_D* -3.0*u + 2.25*u2
  elseif w >=1 and w <=2 then
    var two_minus_u = 2.0 - u
    &w = a_D* 0.25 * (two_minus_u * two_minus_u * two_minus_u)
    &w_dx = a_D* 0.75 * (two_minus_u*two_minus_u)
  else
    &w = 0.0
    &w_dx = 0.0
  end
end

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
  kernel_deval(xi, part1.h, &wi, &wi_dx)
  wi_dx = wi_dx  * inv_hidim_pow_plus_one

  var hj_inv : float = 1.0 / part2.h
  var xj : float= r * hj_inv
  var wj : float = 0.0
  var wj_dx : float = 0.0
  kernel_deval(xj, part2.h, &wj, &wj_dx)
  wj_dx = wj_dx  * inv_hjdim_pow_plus_one

  var r2_eta2 : float = r2 + 0.01*0.01
  var inv_r2eta2 : float = 1.0 / r2_eta2
  var visc_i : float = 4.0 * part1.viscosity
  var visc_j : float = 4.0 * part2.viscosity

  var temp_i : float = visc_i / ( (r2_eta2) * (rhoi + rhoj))
  var temp_j : float = visc_j / ( (r2_eta2) * (rhoi + rhoj))
  var multiplier_i : float = dx0*dx0*wi_dx + dx1*dx1*wi_dx + dx2*dx2*wi_dx
  var multiplier_j : float = dx0*dx0*wj_dx + dx1*dx1*wj_dx + dx2*dx2*wj_dx

  part1.a_hydro_x = part1.a_hydro_x + (mj*temp_i*multiplier_i*dv0)
  part1.a_hydro_y = part1.a_hydro_y + (mj*temp_i*multiplier_i*dv1)
  part1.a_hydro_z = part1.a_hydro_z + (mj*temp_i*multiplier_i*dv2)

  part2.a_hydro_x = part1.a_hydro_x - (mi*temp_j*multiplier_j*dv0)
  part2.a_hydro_y = part1.a_hydro_y - (mi*temp_j*multiplier_j*dv1)
  part2.a_hydro_z = part1.a_hydro_z - (mi*temp_j*multiplier_j*dv2)

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
  
  part1.a_hydro_x = part1.a_hydro_x + hydro0
  part1.a_hydro_y = part1.a_hydro_y + hydro1
  part1.a_hydro_z = part1.a_hydro_z + hydro2

  part2.a_hydro_x = part1.a_hydro_x - hydro0
  part2.a_hydro_y = part1.a_hydro_y - hydro1
  part2.a_hydro_z = part1.a_hydro_z - hydro2

  --Compute velocity gradiants
  var mj_over_rhoj : float = mj * rhoj_inv
  var mi_over_rhoi : float = mi * rhoi_inv
  var dv_over_rho_x_i : float = dv0 * mj_over_rhoj
  var dv_over_rho_y_i : float = dv1 * mj_over_rhoj
  var dv_over_rho_z_i : float = dv2 * mj_over_rhoj
  var dv_over_rho_x_j : float = dv0 * mi_over_rhoi
  var dv_over_rho_y_j : float = dv1 * mi_over_rhoi
  var dv_over_rho_z_j : float = dv2 * mi_over_rhoi

  part1.grad_v_xx = part1.grad_v_xx + (dv_over_rho_x_i * dx0)
  part1.grad_v_xy = part1.grad_v_xy + (dv_over_rho_x_i * dx1)
  part1.grad_v_xz = part1.grad_v_xz + (dv_over_rho_x_i * dx2)
  part1.grad_v_xy = part1.grad_v_xy + (dv_over_rho_y_i * dx0)
  part1.grad_v_yy = part1.grad_v_yy + (dv_over_rho_y_i * dx1)
  part1.grad_v_yz = part1.grad_v_yz + (dv_over_rho_y_i * dx2)
  part1.grad_v_xz = part1.grad_v_xz + (dv_over_rho_z_i * dx0)
  part1.grad_v_yz = part1.grad_v_yz + (dv_over_rho_z_i * dx1)
  part1.grad_v_zz = part1.grad_v_zz + (dv_over_rho_z_i * dx2)

  part2.grad_v_xx = part2.grad_v_xx + (dv_over_rho_x_j * dx0)
  part2.grad_v_xy = part2.grad_v_xy + (dv_over_rho_x_j * dx1)
  part2.grad_v_xz = part2.grad_v_xz + (dv_over_rho_x_j * dx2)
  part2.grad_v_xy = part2.grad_v_xy + (dv_over_rho_y_j * dx0)
  part2.grad_v_yy = part2.grad_v_yy + (dv_over_rho_y_j * dx1)
  part2.grad_v_yz = part2.grad_v_yz + (dv_over_rho_y_j * dx2)
  part2.grad_v_xz = part2.grad_v_xz + (dv_over_rho_z_j * dx0)
  part2.grad_v_yz = part2.grad_v_yz + (dv_over_rho_z_j * dx1)
  part2.grad_v_zz = part2.grad_v_zz + (dv_over_rho_z_j * dx2)

  var acc : float = (part1.pressure + part2.pressure) * ((rhoi_inv * rhoj_inv) * r_inv)

  var dens : float = dv0*wi_dx*dx0 + dv1*wi_dx*dx1 + dv2*wi_dx*dx2

  part1.drho_dt = part1.drho_dt + (mj * dens)
  part1.a_hydro_x = part1.a_hydro_x - (mj * acc * wi_dx * dx0)
  part1.a_hydro_y = part1.a_hydro_y - (mj * acc * wi_dx * dx1)
  part1.a_hydro_z = part1.a_hydro_z - (mj * acc * wi_dx * dx2)

  part2.drho_dt = part2.drho_dt + (mi*dens)
  part2.a_hydro_x = part2.a_hydro_x + (mi * acc * wj_dx * dx0)
  part2.a_hydro_y = part2.a_hydro_y + (mi * acc * wj_dx * dx1)
  part2.a_hydro_z = part2.a_hydro_z + (mi * acc * wj_dx * dx2)

  --Store viscous effect towards CFL condition
  var dvdr : float = dv0*dx0 + dv1*dx1 + dv2*dx2
  var dvdr_rr2 : float = dvdr * inv_r2eta2
  --FIXME : Sort out MAX function
  part1.max_visc = MAX(part1.max_visc, dvdr_rr2)
  part2.max_visc = MAX(part2.max_visc, dvdr_rr2)
  part1.max_visc = MAX(part1.max_visc, dvdr_rr2);

  --FIXME: IMPLEMENT THESE
  [boundary_fluid_interaction( part1, part2, r, r2, dx0, dx1, dx2)];
  [compute_density_diffusive_term(part1, part2, r2, r, wi_dx, wj_dx, dx)];
  [compute_shifting_term(part1, part2, r2, r, wi_dx, wj_dx, dx)];
end
return kernel
end

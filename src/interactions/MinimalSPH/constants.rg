-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

--Cubic spline kernel
local M_PI = 3.141592653589793238462643383279502884
kernel_gamma = global(float, 1.825742)
kernel_constant = global(float, 16.0 / M_PI)
kernel_gamma_inv_dim = global(float, 1.0 / (kernel_gamma:get() * kernel_gamma:get() * kernel_gamma:get()) )
kernel_gamma_inv_dim_plus_one = global(float, 1.0 / (kernel_gamma:get() * kernel_gamma:get() * kernel_gamma:get() * kernel_gamma:get()) )

kernel_const = global(float, 8.0 / 3.0)
--Rough estimate at the moment
--kernel_root = 0.4184291064739227294921875
kernel_root = global(float, 0.5 * kernel_constant:get() * kernel_gamma_inv_dim:get())
--(kernel_coeffs[kernel_degree]) * kernel_constant * \
--   kernel_gamma_inv_dim)
kernel_dimension = global(float, 3.0)
hydro_dimension_inv = global(float, 1.0 / kernel_dimension:get())
hydro_eta = global(float, 1.2348)
hydro_eps = global(float, 1e-4)
alpha = global(float, 0.8)
--Adiabatic index
hydro_gamma = global(float, 5.0/3.0)
hydro_gamma_minus_one = global(float, hydro_gamma:get()-1.0)


--CFL condition
CFL_condition = global(float, 0.1)

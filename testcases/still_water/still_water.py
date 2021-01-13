##-----------------------------------------------------------
## Copyright 2020 Science and Technologies Facilities Council
## Licensed under the MIT License
## Author Aidan Chalk, STFC Hartree Centre

import h5py
import numpy as np

import argparse

pos_x = []
pos_y = []
pos_z = []
cutoff = []
interactions = []
ids = []
count = 0
box_size_x = 25.0
box_size_y = 15.0
box_size_z = 0.0
output_file = 'still_water.hdf5'
WALL=1
FLUID=0

soundspeed_square = 271.247*271.247
gravity = -9.81
abs_grav = 9.81

cs2_rho0_over_7 = soundspeed_square * 1000.0 / 7.0

x_parts = 75
y_parts = 50
#z_parts = 0
particle_spacing = 0.25
box_width = y_parts * particle_spacing
smoothing_len = 1.3 * particle_spacing

parts_per_m = 1.0 / particle_spacing
parts_per_m3 = parts_per_m * parts_per_m * parts_per_m
part_mass = 62.5


diameter = 40
box_bottom = 0.75
box_top = box_bottom + (y_parts) * particle_spacing
box_middle = (box_bottom + box_top) / 2
box_front = 0.5
box_height = box_top - box_bottom
print("height {}".format(box_height))
print("top {}".format(box_top))
box_left = 1.0
box_right = box_left + (x_parts+1) * particle_spacing
print("right {}".format(box_right))

boundary_row = x_parts 
print(boundary_row)
count_boundary = 4 * boundary_row + 8*(2*y_parts+2)

n_hydro = (boundary_row) * (y_parts)

count = count_boundary + n_hydro

boundary_row = boundary_row + 8 #smileyface don't worry its ok )


f = h5py.File(output_file, "w")
dset_posx = f.create_dataset('pos_x', (count,), dtype=np.float64)
dset_posy = f.create_dataset('pos_y', (count,), dtype=np.float64)
dset_cons_acc_y = f.create_dataset('const_acc_y', (count,), dtype=np.float32)
dset_mass = f.create_dataset('mass', (count,), dtype=np.float64)
dset_density = f.create_dataset('density', (count,), dtype=np.float32)
dset_smoothing = f.create_dataset('smoothing_length', (count,), dtype=np.float32)
dset_viscosity = f.create_dataset('viscosity', (count,), dtype=np.float32)
dset_is_bound = f.create_dataset('is_boundary', (count,), dtype=np.int32)

bottom_pressure = 1000.0 * abs_grav * box_height
bottom_density = 1000.0 * pow( 1+(bottom_pressure / cs2_rho0_over_7)   , 1.0/7.0)

boundary_left = box_left - 4*particle_spacing #Starts from 1 index...
##Create bottom
for a in range(boundary_row):
    dset_posx[a] = (boundary_left + particle_spacing * (a+1))
    dset_posy[a] = box_bottom
    dset_cons_acc_y[a] = 0.0
    dset_mass[a] = part_mass
    dset_density[a] = bottom_density
    dset_smoothing[a] = smoothing_len
    dset_viscosity[a] = 1e-6
    dset_is_bound[a] = WALL

    pressure = 1000.0 * abs_grav * (box_height + particle_spacing)
    density = 1000.0 * pow( 1+(pressure / cs2_rho0_over_7)   , 1.0/7.0)
    dset_posx[boundary_row+a] = (boundary_left + particle_spacing * (a+1))
    dset_posy[boundary_row+a] = (box_bottom - particle_spacing)
    dset_cons_acc_y[boundary_row+a] = 0.0
    dset_mass[boundary_row+a] = part_mass
    dset_density[boundary_row+a] = density
    dset_smoothing[boundary_row+a] = smoothing_len
    dset_viscosity[boundary_row+a] = 1e-6
    dset_is_bound[boundary_row+a] = WALL
    
    pressure = 1000.0 * abs_grav * (box_height + 2*particle_spacing)
    density = 1000.0 * pow( 1+(pressure / cs2_rho0_over_7)   , 1.0/7.0)
    dset_posx[2*boundary_row+a] = (boundary_left + particle_spacing * (a+1))
    dset_posy[2*boundary_row+a] = (box_bottom - 2*particle_spacing)
    dset_cons_acc_y[2*boundary_row+a] = 0.0
    dset_mass[2*boundary_row+a] = part_mass
    dset_density[2*boundary_row+a] = density
    dset_smoothing[2*boundary_row+a] = smoothing_len
    dset_viscosity[2*boundary_row+a] = 1e-6
    dset_is_bound[2*boundary_row+a] = WALL
    
    pressure = 1000.0 * abs_grav * (box_height + 3*particle_spacing)
    density = 1000.0 * pow( 1+(pressure / cs2_rho0_over_7)   , 1.0/7.0)
    dset_posx[3*boundary_row+a] = (boundary_left + particle_spacing * (a+1))
    dset_posy[3*boundary_row+a] = (box_bottom - 3*particle_spacing)
    dset_cons_acc_y[3*boundary_row+a] = 0.0
    dset_mass[3*boundary_row+a] = part_mass
    dset_density[3*boundary_row+a] = density
    dset_smoothing[3*boundary_row+a] = smoothing_len
    dset_viscosity[3*boundary_row+a] = 1e-6
    dset_is_bound[3*boundary_row+a] = WALL

#Side Walls
for j in range(1, 2*y_parts):
    pressure = 1000.0 * abs_grav * (box_height - (particle_spacing*(j+1)))
    density = 1000.0 * pow( 1+(pressure / cs2_rho0_over_7)   , 1.0/7.0)
    dset_posx[4*boundary_row + (j-1)*8] = (box_left - particle_spacing*2)
    dset_posy[4*boundary_row + (j-1)*8] = (box_bottom + j*particle_spacing)
    dset_cons_acc_y[4*boundary_row + (j-1)*8] = 0.0
    dset_mass[4*boundary_row + (j-1)*8] = part_mass
    dset_density[4*boundary_row + (j-1)*8] = density
    dset_smoothing[4*boundary_row + (j-1)*8] = smoothing_len
    dset_viscosity[4*boundary_row + (j-1)*8] = 1e-6
    dset_is_bound[4*boundary_row + (j-1)*8] = WALL

    dset_posx[4*boundary_row + (j-1)*8 +1] = (box_left - particle_spacing)
    dset_posy[4*boundary_row + (j-1)*8+1] = (box_bottom + j*particle_spacing)
    dset_cons_acc_y[4*boundary_row + (j-1)*8+1] = 0.0
    dset_mass[4*boundary_row + (j-1)*8+1] = part_mass
    dset_density[4*boundary_row + (j-1)*8+1] = density
    dset_smoothing[4*boundary_row + (j-1)*8+1] = smoothing_len
    dset_viscosity[4*boundary_row + (j-1)*8+1] = 1e-6
    dset_is_bound[4*boundary_row + (j-1)*8+1] = WALL
    
    dset_posx[4*boundary_row + (j-1)*8 +2] = (box_left - 3*particle_spacing)
    dset_posy[4*boundary_row + (j-1)*8 +2] = (box_bottom + j*particle_spacing)
    dset_cons_acc_y[4*boundary_row + (j-1)*8 +2] = 0.0
    dset_mass[4*boundary_row + (j-1)*8 +2] = part_mass
    dset_density[4*boundary_row + (j-1)*8 +2] = density
    dset_smoothing[4*boundary_row + (j-1)*8 +2] = smoothing_len
    dset_viscosity[4*boundary_row + (j-1)*8 +2] = 1e-6
    dset_is_bound[4*boundary_row + (j-1)*8 +2] = WALL
    
    dset_posx[4*boundary_row + (j-1)*8 +3] = (box_left)
    dset_posy[4*boundary_row + (j-1)*8 +3] = (box_bottom + j*particle_spacing)
    dset_cons_acc_y[4*boundary_row + (j-1)*8 +3] = 0.0
    dset_mass[4*boundary_row + (j-1)*8 +3] = part_mass
    dset_density[4*boundary_row + (j-1)*8 +3] = density
    dset_smoothing[4*boundary_row + (j-1)*8 +3] = smoothing_len
    dset_viscosity[4*boundary_row + (j-1)*8 +3] = 1e-6
    dset_is_bound[4*boundary_row + (j-1)*8 +3] = WALL
    
    dset_posx[4*boundary_row + (j-1)*8 +4] = (box_right)
    dset_posy[4*boundary_row + (j-1)*8 +4] = (box_bottom + j*particle_spacing)
    dset_cons_acc_y[4*boundary_row + (j-1)*8 +4] = 0.0
    dset_mass[4*boundary_row + (j-1)*8 +4] = part_mass
    dset_density[4*boundary_row + (j-1)*8 +4] = density
    dset_smoothing[4*boundary_row + (j-1)*8 +4] = smoothing_len
    dset_viscosity[4*boundary_row + (j-1)*8 +4] = 1e-6
    dset_is_bound[4*boundary_row + (j-1)*8 +4] = WALL
    
    dset_posx[4*boundary_row + (j-1)*8 +5] = (box_right + particle_spacing)
    dset_posy[4*boundary_row + (j-1)*8 +5] = (box_bottom + j*particle_spacing)
    dset_cons_acc_y[4*boundary_row + (j-1)*8 +5] = 0.0
    dset_mass[4*boundary_row + (j-1)*8 +5] = part_mass
    dset_density[4*boundary_row + (j-1)*8 +5] = density
    dset_smoothing[4*boundary_row + (j-1)*8 +5] = smoothing_len
    dset_viscosity[4*boundary_row + (j-1)*8 +5] = 1e-6
    dset_is_bound[4*boundary_row + (j-1)*8 +5] = WALL
    
    dset_posx[4*boundary_row + (j-1)*8 +6] = (box_right + 2*particle_spacing)
    dset_posy[4*boundary_row + (j-1)*8 +6] = (box_bottom + j*particle_spacing)
    dset_cons_acc_y[4*boundary_row + (j-1)*8 +6] = 0.0
    dset_mass[4*boundary_row + (j-1)*8 +6] = part_mass
    dset_density[4*boundary_row + (j-1)*8 +6] = density
    dset_smoothing[4*boundary_row + (j-1)*8 +6] = smoothing_len
    dset_viscosity[4*boundary_row + (j-1)*8 +6] = 1e-6
    dset_is_bound[4*boundary_row + (j-1)*8 +6] = WALL
    
    dset_posx[4*boundary_row + (j-1)*8 +7] = (box_right + 3*particle_spacing)
    dset_posy[4*boundary_row + (j-1)*8 +7] = (box_bottom + j*particle_spacing)
    dset_cons_acc_y[4*boundary_row + (j-1)*8 +7] = 0.0
    dset_mass[4*boundary_row + (j-1)*8 +7] = part_mass
    dset_density[4*boundary_row + (j-1)*8 +7] = density
    dset_smoothing[4*boundary_row + (j-1)*8 +7] = smoothing_len
    dset_viscosity[4*boundary_row + (j-1)*8 +7] = 1e-6
    dset_is_bound[4*boundary_row + (j-1)*8 +7] = WALL

##Boundaries done I hope 
for i in range(0, x_parts):
    for j in range(0, y_parts):
        pressure = 1000.0 * abs_grav * (box_height - (particle_spacing*(j+1)))
        density = 1000.0 * pow( 1+(pressure / cs2_rho0_over_7)   , 1.0/7.0)
        dset_posx[count_boundary + i*(y_parts)+j] = (box_left + particle_spacing*(i+1))
        dset_posy[count_boundary + i*(y_parts)+j] = (box_bottom + particle_spacing*(j+1))
        dset_cons_acc_y[count_boundary + i*(y_parts)+j] = gravity
        dset_mass[count_boundary + i*(y_parts)+j] = part_mass
        dset_density[count_boundary + i*(y_parts)+j] = density
        dset_smoothing[count_boundary + i*(y_parts)+j] = smoothing_len
        dset_viscosity[count_boundary + i*(y_parts)+j] = 1e-6
        dset_is_bound[count_boundary + i*(y_parts)+j] = FLUID

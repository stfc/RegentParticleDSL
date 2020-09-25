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
x_count = 0
y_count = 0
z_count = 0
box_size_x = 3.0
box_size_y = 3.0
box_size_z = 3.0
set_cutoff = 1.75
output_file = ''

def parse_inputs():
    global box_size_x, box_size_y, box_size_z, x_count, y_count, z_count, output_file, set_cutoff
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, type=str, help='Path of the output file')
    parser.add_argument('--x_count', required=True, type=int, help='Number of particles in x dimension')
    parser.add_argument('--y_count', required=True, type=int, help='Number of particles in y dimension')
    parser.add_argument('--z_count', required=True, type=int, help='Number of particles in z dimension')
    parser.add_argument('--cutoff', required=False, type=float, help='Size of the cutoff radius for the particles in the simulation')
    parser = parser.parse_args()
    output_file = parser.output
    x_count = parser.x_count
    y_count = parser.y_count
    z_count = parser.z_count
    box_size_x = float(x_count)
    box_size_y = float(y_count)
    box_size_z = float(z_count)
    if(parser.cutoff):
        set_cutoff = parser.cutoff

def create_problem():
    global x_count, y_count, z_count, output_file, set_cutoff
    f = h5py.File(output_file, "w")
    count = x_count * y_count * z_count
    dset_posx = f.create_dataset('pos_x', (count,), dtype=np.float64)
    dset_posy = f.create_dataset('pos_y', (count,), dtype=np.float64)
    dset_posz = f.create_dataset('pos_z', (count,), dtype=np.float64)
    dset_cutoff = f.create_dataset('cutoff', (count,), dtype=np.float64)
    dset_interactions = f.create_dataset('interactions', (count,), dtype=np.uint32)
    dset_ids = f.create_dataset('ids', (count,), dtype=np.int64)
    index = 0
    for x in range(x_count):
        for y in range(y_count):
            for z in range(z_count):
                dset_posx[index] = x
                dset_posy[index] = y
                dset_posz[index] = z
                index = index + 1
    for i in range(count):
        dset_cutoff[i] = set_cutoff
        dset_interactions[i] = 0
        dset_ids[i] = i

def main():
    parse_inputs()
    create_problem()

main()
        



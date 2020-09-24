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
box_size_x = 3.0
box_size_y = 3.0
box_size_z = 3.0
set_cutoff = 1.75
output_file = ''
seed = 1

def parse_inputs():
    global count, output_file, set_cutoff, seed
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, type=str, help='Path of the output file')
    parser.add_argument('--count', required=True, type=int, help='Number of particles')
    parser.add_argument('--seed', required=True, type=int, help='Random seed for the test')
    parser = parser.parse_args()
    output_file = parser.output
    count = parser.count
    seed = parser.seed

def create_problem():
    global count, output_file, set_cutoff
    f = h5py.File(output_file, "w")
    dset_posx = f.create_dataset('pos_x', (count,), dtype=np.float64)
    dset_posy = f.create_dataset('pos_y', (count,), dtype=np.float64)
    dset_posz = f.create_dataset('pos_z', (count,), dtype=np.float64)
    dset_cutoff = f.create_dataset('cutoff', (count,), dtype=np.float64)
    dset_interactions = f.create_dataset('interactions', (count,), dtype=np.uint32)
    dset_ids = f.create_dataset('ids', (count,), dtype=np.int64)
    index = 0
    cutoff = np.random.rand(1)
    x_vals = np.random.rand(count) * 3.0
    y_vals = np.random.rand(count) * 3.0
    z_vals = np.random.rand(count) * 3.0
    set_cutoff = cutoff[0]
    for x in range(count):
        dset_posx[x] = x_vals[x]
        dset_posy[x] = y_vals[x]
        dset_posz[x] = z_vals[x]
    for i in range(count):
        dset_cutoff[i] = set_cutoff
        dset_interactions[i] = 0
        dset_ids[i] = i

def main():
    parse_inputs()
    create_problem()

main()
        



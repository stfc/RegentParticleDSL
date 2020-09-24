import h5py
import numpy as np

import argparse

pos_x = []
pos_y = []
pos_z = []
cutoff = []
interactions = []
ids = []
box_size_x = 3.0
box_size_y = 3.0
box_size_z = 3.0
input_file = ''
output_file = ''

def parse_inputs():
    global box_size_x, box_size_y, box_size_z, input_file, output_file
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, type=str, help='Path of the input file')
    parser.add_argument('--output', required=True, type=str, help='Path of the output file')
    parser.add_argument('--box_x', required=True, type=float)
    parser.add_argument('--box_y', required=True, type=float)
    parser.add_argument('--box_z', required=True, type=float)
    parser = parser.parse_args()
    input_file = parser.input
    output_file = parser.output
    box_size_x = parser.box_x
    box_size_y = parser.box_y
    box_size_z = parser.box_z


def read_file(filename):
    global pos_x, pos_y, pos_z, cutoff, interactions, ids
    f = h5py.File(filename, 'r')
    print(list(f.keys()))
    print(f['pos_x'].dtype)
    print(f['cutoff'].dtype)
    print(f['interactions'].dtype)
    print(f['ids'].dtype)
    pos_x = np.zeros(f['pos_x'].shape, dtype = f['pos_x'].dtype)
    pos_y = np.zeros(f['pos_y'].shape, dtype = f['pos_y'].dtype)
    pos_z = np.zeros(f['pos_z'].shape, dtype = f['pos_z'].dtype)
    cutoff = np.zeros(f['cutoff'].shape, dtype = f['cutoff'].dtype)
    interactions = np.zeros(f['interactions'].shape, dtype = f['interactions'].dtype)
    ids = np.zeros(f['ids'].shape, dtype = f['ids'].dtype)
    for i in range(len(pos_x)):
        pos_x[i] = f['pos_x'][i]
        pos_y[i] = f['pos_y'][i]
        pos_z[i] = f['pos_z'][i]
        cutoff[i] = f['cutoff'][i]
        ids[i] = f['ids'][i]

def write_file(filename):
    global pos_x, pos_y, pos_z, cutoff, interactions, ids
    f = h5py.File(filename, "w")
    dset_posx = f.create_dataset('pos_x', pos_x.shape, dtype=pos_x.dtype)
    dset_posy = f.create_dataset('pos_y', pos_y.shape, dtype=pos_y.dtype)
    dset_posz = f.create_dataset('pos_z', pos_z.shape, dtype=pos_z.dtype)
    dset_cutoff = f.create_dataset('cutoff', cutoff.shape, dtype=cutoff.dtype)
    dset_interactions = f.create_dataset('interactions', interactions.shape, dtype=interactions.dtype)
    dset_ids = f.create_dataset('ids', ids.shape, dtype=ids.dtype)
    for i in range(len(pos_x)):
        dset_posx[i] = pos_x[i]
        dset_posy[i] = pos_y[i]
        dset_posz[i] = pos_z[i]
        dset_cutoff[i] = cutoff[i]
        dset_interactions[i] = interactions[i]
        dset_ids[i] = ids[i]
    print(interactions)
    print(filename)

def n2_interactions():
    global pos_x, pos_y, pos_z, cutoff, interactions, box_size_x, box_size_y, box_size_z
    half_box_x = 0.5 * box_size_x
    half_box_y = 0.5 * box_size_y
    half_box_z = 0.5 * box_size_z
    for i in range(len(pos_x)):
        for j in range(len(pos_x)):
            if i >= j:
                continue
            dx = pos_x[i] - pos_x[j]
            dy = pos_y[i] - pos_y[j]
            dz = pos_z[i] - pos_z[j]
            if dx > half_box_x:
                dx = dx - box_size_x
            if dy > half_box_y:
                dy = dy - box_size_y
            if dz > half_box_z:
                dz = dz - box_size_z
            if dx < -half_box_x: 
                dx = dx + box_size_x
            if dy < -half_box_y: 
                dy = dy + box_size_y
            if dz < -half_box_z: 
                dz = dz + box_size_z
            cutoff_i = max(cutoff[i], cutoff[j])
            r2 = dx*dx + dy*dy + dz*dz
            cutoff2 = cutoff_i*cutoff_i
            if r2 <= cutoff2:
                interactions[i] = interactions[i] + 1
                interactions[j] = interactions[j] + 1

def main():
    global input_file, output_file
    parse_inputs()
    read_file(input_file)
    n2_interactions()
    write_file(output_file)

main()

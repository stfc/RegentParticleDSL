##-----------------------------------------------------------
## Copyright 2020 Science and Technologies Facilities Council
## Licensed under the MIT License
## Author Aidan Chalk, STFC Hartree Centre

import h5py
import numpy as np
import matplotlib.pyplot as plt
import argparse




#f = h5py.File('still_water.hdf5', 'r')
#f = h5py.File('initial_output.hdf5', 'r')
f = h5py.File('outputs/file9.hdf5', 'r')

#plt.scatter(f['pos_x'], f['pos_y'], c=f['is_boundary'])
print(max(f['pos_x']))
print(max(f['pos_y']))
plt.scatter(f['pos_x'], f['pos_y'], c=f['density'])
plt.show()

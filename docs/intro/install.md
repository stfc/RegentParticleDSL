# Installing HartreeParticleDSL

## Downloading the repository
The repository is available on [github here](https://github.com/stfc/HartreeParticleDSL). The master branch contains the current version of the code, however this can be modified regularly, so if there are issues please let us know and try the tagged releases instead.

## Requirements

HartreeParticleDSL requires nothing special other than a working version of Regent, and a HDF5 installation.
The easiest way to have this is to use the Dockerfile available in the repository, and mount the DSL into the docker image. 
On linux, this will build and run the docker image, and mount the currently directory as `/DSL`. The final line shows how to launch a
code implementation:
```
  docker build -t regentparticle
  docker run -it -v $PWD:/DSL regentparticle
  cd DSL
  regent /path/to/code.rg -fflow 0
```

Alternatively, these libraries can be built yourself on linux. Instructions for Regent are available on [regent-lang.org/install](http://regent-lang.org/install/), 
and HDF5 is available at [hdf5group.org](https://www.hdfgroup.org/solutions/hdf5/). Regent needs to be built with the `--hdf5` option. To build Regent with HDF5,
make sure the appropriate paths to HDF5 are visible on the `PATH`, `LD_LIBRARY_PATH`, `LIBRARY_PATH`, `INCLUDE_PATH` and `CPLUS_INCLUDE_PATH`.


## Project Outline
The main HartreeParticleDSL folder contains 3 subdirectories. The `example_types` directory contains some outdated example code, which 
has been superseded at the moment by this documentation.

The `examples` folder contains small test programs, that have working programs that can be used to check the installation has worked correctly.

The `src` folder contains all the DSL library code, such as neighbour search algorithms, default type declarations and some kernel implementations.

At the moment, all HartreeParticleDSL programs should be launched from the root directory of the DSL, due to how files are imported.

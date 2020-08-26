#RegentParticleDSL
#This includes HDF5

FROM ubuntu:14.04

MAINTAINER Aidan Chalk <aidan.chalk@stfc.ac.uk>

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && \
    apt-get install -y build-essential clang-3.5 git libclang-3.5-dev libncurses5-dev llvm-3.5-dev wget zlib1g-dev libhdf5-dev cmake3 && \
    apt-get clean
    
RUN git clone https://github.com/StanfordLegion/legion.git

RUN cd legion/language
RUN LLVM_CONFIG=llvm-config-3.5 python3 legion/language/install.py --rdir=auto --hdf5 --no-terra-cmake
RUN ln -s /legion/language/regent.py /usr/local/bin/regent
    
CMD ["/bin/bash"]

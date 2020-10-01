#RegentParticleDSL
#This includes HDF5

FROM ubuntu:14.04

MAINTAINER Aidan Chalk <aidan.chalk@stfc.ac.uk>

ARG BRANCH_NAME
ARG COMMIT_HASH
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && \
    apt-get install -y build-essential clang-3.5 git libclang-3.5-dev libncurses5-dev llvm-3.5-dev wget zlib1g-dev libhdf5-dev cmake3 vim && \
    apt-get clean
    
RUN git clone -b master https://github.com/StanfordLegion/legion.git

RUN cd legion/language
RUN LLVM_CONFIG=llvm-config-3.5 python3 legion/language/install.py --rdir=auto --hdf5 --no-terra-cmake
RUN ln -s /legion/language/regent.py /usr/local/bin/regent
RUN cd /
RUN git clone https://github.com/StanfordLegion/regent.vim.git
RUN mkdir ~/.vim/
RUN mkdir ~/.vim/syntax/
RUN cp regent.vim/regent.vim ~/.vim/syntax/.
RUN echo "au BufNewFile,BufRead *.rg set filetype=regent" >> ~/.vimrc
RUN rm -rf regent.vim
RUN git clone -b $BRANCH_NAME https://github.com/stfc/RegentParticleDSL.git
RUN cd RegentParticleDSL && git checkout $COMMIT_HASH

RUN wget https://www.python.org/ftp/python/3.6.3/Python-3.6.3.tgz
RUN sudo apt install -y libssl-dev libncurses5-dev libsqlite3-dev libreadline-dev libtk8.5 libgdm-dev libdb4o-cil-dev libpcap-dev
RUN tar -xvf Python-3.6.3.tgz
RUN cd Python-3.6.3 && sudo ./configure && sudo make -j2 && sudo make install
RUN python3.6 -m pip install numpy h5py

ENV HDF5_INCLUDE_PATH=/usr/include/
    
CMD ["/bin/bash"]

#!/bin/bash

set -euo pipefail

cd libsbml-5.20.5 || exit 1
./configure --enable-fbc 

if [ ! -d build ]; then
  mkdir build
fi
cd build || exit 1
mkdir libsbml 
installpath=$PWD/libsbml

cmake ..\
    -DCMAKE_INSTALL_PREFIX=$installpath \
    -DENABLE_FBC=ON 

make -j4
make install

# TODO: fix this if possible
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH$installpath/lib">> ~/.bashrc
source ~/.bashrc
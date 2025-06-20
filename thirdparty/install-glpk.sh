#!/bin/bash

set -euo pipefail

cd glpk-4.35 || exit 1

if [ ! -d glpk ]; then
  mkdir glpk
fi
installpath=$PWD/glpk

./configure --prefix=$installpath

make
make install

# TODO: fix this if possible
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$installpath/lib:">> ~/.bashrc
source ~/.bashrc
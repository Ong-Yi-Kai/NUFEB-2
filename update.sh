#!/bin/bash

set -euo pipefail

cd ${0%/*} || exit 1 # Run from this directory

BASEDIR="./src/NUFEB"
LAMMPS_DIR="./lammps_stable_23Jun2022/src"

# read the filename from command line to be update
if [ $# -ne 1 ]; then
    echo "Requires filename. Usage: '$0 <filename>'"
    exit 1
fi

filepath_h="$BASEDIR/$1.h"
filepath_cpp="$BASEDIR/$1.cpp"

# check if the file exists
if [ ! -f "$filepath_h" ]; then
    echo "File $filepath_h does not exist."
    exit 1
fi
if [ ! -f "$filepath_cpp" ]; then
    echo "File $filepath_cpp does not exist."
    exit 1
fi

# copy the .cpp and .h file over to LAMMPS source
new_filepath_h="$LAMMPS_DIR/NUFEB/$1.h"
cp $filepath_h $new_filepath_h
new_filepath_cpp="$LAMMPS_DIR/NUFEB/$1.cpp"
cp $filepath_cpp $new_filepath_cpp

# make the build system happy
cd $LAMMPS_DIR || exit 1
make yes-nufeb

# build the new executable
make -j4 mpi

# move the new executable to the current directory
cd ../..
mv ./lammps_stable_23Jun2022/src/lmp_mpi ./nufeb_mpi
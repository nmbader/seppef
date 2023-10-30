#!/usr/bin/env bash
# inspired by build script for Arch Linux fftw pacakge:
mkdir build
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX
make
make install

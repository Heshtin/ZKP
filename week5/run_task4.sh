#!/bin/bash

set -e

cd /mnt/c/Users/bluni/ZKP

export LD_LIBRARY_PATH=/mnt/c/Users/bluni/ZKP/mcl/lib:$LD_LIBRARY_PATH

echo "Compiling Task4.cpp with Task3.cpp..."

g++ Task4.cpp Task3.cpp \
    -I ./mcl/include \
    -L ./mcl/lib \
    -lmcl \
    -std=c++17 \
    -O2 \
    -o task4

echo "Compilation successful."
echo "Running Task 4..."
echo

./task4
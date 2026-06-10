#!/bin/bash

set -e

cd /mnt/c/Users/bluni/ZKP

export LD_LIBRARY_PATH=/mnt/c/Users/bluni/ZKP/mcl/lib:$LD_LIBRARY_PATH

echo "Compiling Task3.cpp..."

g++ Task3.cpp \
    -I ./mcl/include \
    -L ./mcl/lib \
    -lmcl \
    -std=c++17 \
    -O2 \
    -o task3

echo "Compilation successful."
echo "Running Task 3..."
echo

./task3
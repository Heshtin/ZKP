#!/bin/bash

set -e

echo "Compiling R1CS example..."

g++ R1CS.cpp \
    -I ./mcl/include \
    -L ./mcl/lib \
    -lmcl \
    -std=c++17 \
    -O2 \
    -o r1cs

echo "Compilation successful."
echo

LD_LIBRARY_PATH=./mcl/lib ./r1cs
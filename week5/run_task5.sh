#!/bin/bash

set -e

echo "Compiling Task 5: Univariate SumCheck..."

g++ Task5.cpp Task3.cpp \
    -I ./mcl/include \
    -L ./mcl/lib \
    -lmcl \
    -std=c++17 \
    -O2 \
    -o task5

echo "Compilation successful."
echo "Running Task 5..."
echo

LD_LIBRARY_PATH=./mcl/lib ./task5
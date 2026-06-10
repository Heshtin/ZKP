#!/bin/bash
echo "Compiling Task 1: NTT & INTT"
g++ Task1_mcl.cpp \
    -I ./mcl/include \
    -L ./mcl/lib \
    -lmcl \
    -std=c++17 \
    -O2 \
    -o task1_mcl

export LD_LIBRARY_PATH=./mcl/lib:$LD_LIBRARY_PATH
echo "Compilation successful."
echo "Running Task 1..."
echo

./task1_mcl
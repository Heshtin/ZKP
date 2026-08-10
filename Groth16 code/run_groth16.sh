#!/bin/bash

set -e

echo "Compiling Groth16 setup "

mkdir -p output

g++ main_groth16.cpp \
    Groth16.cpp \
    QAP.cpp \
    R1CS.cpp \
    -I ./mcl/include \
    -L ./mcl/lib \
    -lmcl \
    -std=c++17 \
    -O2 \
    -o output/groth16_app

echo "Compilation successful."
echo "Running Groth16 setup"
echo

LD_LIBRARY_PATH=./mcl/lib ./output/groth16_app
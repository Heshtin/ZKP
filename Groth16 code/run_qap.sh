#!/bin/bash

set -e

echo "Compiling QAP example..."

g++ main_qap.cpp QAP.cpp R1CS.cpp \
    -I ./mcl/include \
    -L ./mcl/lib \
    -lmcl \
    -std=c++17 \
    -O2 \
    -o qap

echo "Compilation successful."
echo

LD_LIBRARY_PATH=./mcl/lib ./qap
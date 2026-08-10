#pragma once

#include <vector>
#include <mcl/bn.hpp>

using std::vector;
using mcl::bn::Fr;

struct R1CS {
    vector<vector<Fr>> A;
    vector<vector<Fr>> B;
    vector<vector<Fr>> C;
    int v;
    int PI;
};


Fr multiply(
    const vector<Fr>& row,
    const vector<Fr>& assign
);
bool verify_R1CS (
    const R1CS& r1cs,
    const vector<Fr>& assign
);

R1CS poly_R1CS (
    vector<Fr> coeff
);

vector<Fr> witness_vector(
    vector<Fr> coeff,
    const Fr& x
);
Fr random_fr();

vector<Fr> rand_poly(int deg);

void print_poly(
    const vector<Fr>& coeff
);
#pragma once

#include "R1CS.hpp"

#include <vector>
#include <utility>

using std::pair;
using std::vector;
using mcl::bn::Fr;

struct QAP {
    vector<vector<Fr>> u;
    vector<vector<Fr>> v;
    vector<vector<Fr>> w;

    vector<Fr> t;
    vector<Fr> pt;

    int num_v;
    int num_c;
    int PI;
};

Fr eval_poly(
    const vector<Fr>& poly,
    const Fr& x
);

QAP convert(
    const R1CS& r1cs
);

bool verify_QAP(
    const QAP& qap,
    const vector<Fr>& assig,
    vector<Fr>& h
);

void print_checks(
    const QAP& qap,
    const vector<Fr>& assign
);
#pragma once

#include "QAP.hpp"

#include <vector>
#include <mcl/bn.hpp>

using std::vector;

using mcl::bn::Fr;
using mcl::bn::G1;
using mcl::bn::G2;

struct PKey {
    G1 a_G1;
    G1 b_G1;
    G2 b_G2;
    G1 d_G1;
    G2 d_G2;

    vector<G1> uQ_G1;
    vector<G1> vQ_G1;
    vector<G2> vQ_G2;

    vector<G1> private_G1;
    vector<G1> hQ_G1;

    int num_v;
    int PI;
};

struct VKey {
    G1 a_G1;
    G2 b_G2;
    G2 g_G2;
    G2 d_G2;

    vector<G1> public_G1;

    int PI;
};

struct Keys {
    PKey pKey;
    VKey vKey;
};

struct GrothProof {
    G1 A;
    G2 B;
    G1 C;
};

Fr rand_fr();

Keys groth_Setup(
    const QAP& qap
);

GrothProof groth_Prove(
    const PKey& pKey,
    const QAP& qap,
    const vector<Fr>& assign
);

bool groth_Verify(
    const VKey& vKey,
    const vector<Fr>& PI,
    const GrothProof& proof
);
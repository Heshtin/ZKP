#ifndef KZG_HPP
#define KZG_HPP

#include <iostream>
#include <vector>
#include <mcl/bn.hpp>

using namespace std;
using namespace mcl::bn;

struct PK {
    int t;
    vector<G1> gp;
    G2 h;
    G2 tauH;
};

struct Witness {
    Fr i;
    Fr value;
    G1 w;
};
void init();

Fr random_fr();

PK Setup(int maxDegree);
G1 Commit(const vector<Fr>& coeffs, const PK& pk);

vector<Fr> Open(const PK& pk, const G1& C, const vector<Fr>& poly);
bool VerifyPoly(const PK& pk, const G1& C, const vector<Fr>& poly);

Fr evalPoly(const vector<Fr>& f, const Fr& x);

vector<Fr> divideByLinear(vector<Fr> f, const Fr& point, const Fr& y);

Witness CreateWitness(const vector<Fr>& f, const Fr& i, const PK& pk);
bool VerifyEval(const PK& pk, const G1& C, const Witness& W);

#endif
#include "Task3.hpp"

struct SumCheckProof {
    bool valid = false;
    G1 Cp;
    G1 Ch;
    G1 Cf;

    Witness Wp;
    Witness Wh;
    Witness Wf;

    //vector<Fr> hPoly;
    //vector<Fr> fPoly;
};

void trim(vector<Fr>& poly) {
    while (poly.size() > 1 && poly.back().isZero()) {
        poly.pop_back();
    }
}

pair<vector<Fr>, vector<Fr>> divide_vanish_Poly (vector<Fr> p, int n) {
    trim(p);

    int degree = (int)p.size() - 1;
    if (degree < n) {
        return {{0}, p};
    }

    vector<Fr> h(degree - n + 1);
    vector<Fr> r = p;
    for (int k = degree; k >= n; k--) {
        Fr temp = r[k];
        if (!temp.isZero()) {
            h[k - n] = temp;
            r[k] -= temp;
            r[k - n] += temp;
        }
    }

    r.resize(n);
    trim(h);
    trim(r);

    return {h, r};
}
vector<Fr> divideByX(const vector<Fr>& r) {
    if (r.size() <= 1) {
        return {0};
    }
    vector<Fr> f(r.size() - 1);
    for (int i = 1; i < (int)r.size(); i++) {
        f[i - 1] = r[i];
    }
    trim(f);
    return f;
}
//Prover
SumCheckProof Prover_proof(const vector<Fr>& p,int n, const Fr& r,const PK& pk) {
    SumCheckProof proof;

    proof.Cp = Commit(p, pk);

    //p = h * Z_H + r
    auto [h, rem] = divide_vanish_Poly(p, n);

    if (!rem[0].isZero()) {
        cout <<"SumCheck failed: remainder constant term is not zero."<<endl;
        return proof;
    }
    proof.valid = true;

    // r = X * f
    vector<Fr> f = divideByX(rem);

    //proof.hPoly = h;
    //proof.fPoly = f;

    proof.Ch = Commit(h, pk);
    proof.Cf = Commit(f, pk);

    proof.Wp = CreateWitness(p, r, pk);
    proof.Wh = CreateWitness(h, r, pk);
    proof.Wf = CreateWitness(f, r, pk);

    return proof;
}
//Verifier
bool Verifier_SumCheck(
    const SumCheckProof& proof,
    int n,
    const Fr& r,
    const PK& pk) {
    if (!proof.valid) {
        return false;
    }
    bool okP = VerifyEval(pk, proof.Cp, proof.Wp);
    bool okH = VerifyEval(pk, proof.Ch, proof.Wh);
    bool okF = VerifyEval(pk, proof.Cf, proof.Wf);

    if (!okP || !okH || !okF) {
        return false;
    }

    // p(r) = h(r) *(r^n-1)+ r*f(r)
    Fr ZH;
    Fr::pow(ZH , r, n);
    ZH -= 1;
    Fr rhs = proof.Wh.value;
    rhs *= ZH;
    Fr rF = r;
    rF *= proof.Wf.value;
    rhs += rF;
    return proof.Wp.value == rhs;
}
void testSumCheck(const vector<Fr>& p, int n) {
    PK pk = Setup((int)p.size() - 1);
    Fr r = random_fr();

    SumCheckProof proof = Prover_proof(p, n, r, pk);

    if (!proof.valid) {
        cout << "Univariate SumCheck test failed." << endl;
        return;
    }
    bool ok = Verifier_SumCheck(proof, n, r, pk);

    cout << "Univariate SumCheck test" << endl;
    cout << "n = " << n << endl;
    cout << "r = " << r << endl;
    cout << "p(r) = " << proof.Wp.value << endl;
    cout << "h(r) = " << proof.Wh.value << endl;
    cout << "f(r) = " << proof.Wf.value << endl;
    cout << "VerifySumCheck: " << (ok ? "passed" : "failed") << endl;
}
int main() {
    init();
    // h(X) = 2 + 3X
    // f(X) = 5 + 6X + 7X^2
    //
    // p(X) = h(X)(X^4 - 1) + Xf(X)

    vector<Fr> h = {2, 3};
    vector<Fr> f = {5, 6, 7};

    vector<Fr> p ={-2, 2, 6, 7, 2, 3};

    testSumCheck(p,4);
    return 0;
}
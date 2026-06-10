#include "Task3.hpp"

void init() {
    initPairing(BN_SNARK1);
}
// random element.
Fr random_fr() {
    Fr x;
    x.setByCSPRNG();
    return x;
}
PK Setup(int maxDegree) {
    PK pk;
    pk.t = maxDegree;
    G1 g;
    G2 h;

    mapToG1(g, 1);
    mapToG2(h, 1);

    pk.h = h;
    Fr tau = random_fr();
    Fr power = 1;
    for (int i = 0; i <= maxDegree; i++) {
        G1 temp;
        G1::mul(temp, g, power);
        pk.gp.push_back(temp);

        power *= tau;
    }
    G2::mul(pk.tauH, h, tau);
    return pk;
}
G1 Commit(const vector<Fr>& coeffs,const PK& pk) {
    G1 C;
    C.clear();
    for (int i = 0; i < coeffs.size(); i++) {
        G1 temp;
        G1::mul(temp, pk.gp[i], coeffs[i]);
        G1::add(C, C, temp);
    }
    return C;
}
vector<Fr> Open(const PK& pk,const G1& C,const vector<Fr>& poly) {
    return poly;
}
bool VerifyPoly(const PK& pk,const G1& C, const vector<Fr>& poly) {
    if ((int)poly.size()> pk.t + 1) {
        return false;
    }
    G1 resC;
    resC.clear();
    for (int i = 0; i < (int)poly.size(); i++) {
        G1 temp;
        G1::mul(temp, pk.gp[i], poly[i]);
        G1::add(resC, resC, temp);
    }
    return resC == C;
}
Fr evalPoly(const vector<Fr>& f, const Fr& x) {
    Fr result = 0;
    for (int i = f.size() - 1; i >= 0; i--) {
        result *= x;
        result += f[i];
    }
    return result;
}
vector<Fr> divideByLinear(vector<Fr> f, const Fr& i, const Fr& y) {
    int n = f.size();
    if (n == 1) {
        return {0};
    }
    f[0] -= y;
    vector<Fr> q(n - 1);

    q[n - 2] = f[n - 1];
    for (int j = n - 3; j >= 0; j--) {
        q[j] = f[j + 1] + i * q[j + 1];
    }

    return q;
}

Witness CreateWitness(const vector<Fr>& f, const Fr& i, const PK& pk) {
    Witness witness;
    witness.i = i;

    Fr y = evalPoly(f, i);
    witness.value = y;

    vector<Fr> q = divideByLinear(f, i, y);
    witness.w = Commit(q, pk);

    return witness;
}
bool VerifyEval(const PK& pk, const G1& C, const Witness& W) {
    G1 yG;
    G1::mul(yG, pk.gp[0], W.value);

    G1 leftG1;
    G1::sub(leftG1, C, yG);

    G2 iH;
    G2::mul(iH, pk.h, W.i);

    G2 rightG2;
    G2::sub(rightG2, pk.tauH, iH);

    Fp12 e1, e2;

    pairing(e1, leftG1, pk.h);
    pairing(e2, W.w, rightG2);

    return e1 == e2;
}
void test(const vector<Fr>& f, const Fr& i){
    PK pk = Setup((int)f.size() - 1);
    G1 C = Commit(f, pk);

    vector<Fr> openedPoly = Open(pk, C, f);
    bool polyOk = VerifyPoly(pk, C, openedPoly);

    Witness W = CreateWitness(f, i, pk);
    bool ok = VerifyEval(pk, C, W);
    cout << "Example: \n" << endl;
    cout << "f(i) = " << W.value << endl;
    cout << "i = " << i <<"\nsize = " << pk.t << endl;
    cout << "VerifyPoly: " << (polyOk ? "passed" : "failed") << endl;
    cout << "VerifyEval: " << (ok ? "passed" : "failed") << endl;
    cout << endl;
}
// int main() {
//     init();

//     // f(x) = 3 + 2x + x^2
//     vector<Fr> f1 = {3,2,1};
//     vector<Fr> f2 = {10,9,8,4,3,2,1};
    
//     test(f1,4);
//     test(f1,random_fr());
//     test(f2,random_fr());

//     return 0;
// }
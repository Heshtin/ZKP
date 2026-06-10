#include "Task3.hpp"
struct ZeroTestProof {
    bool valid = false;
    PK pk;
    G1 Cf;       
    G1 Cq;       
    Witness Wf;
    Witness Wq;
};
vector<Fr> vanishPoly(vector<Fr> f, int size) {
    int n = f.size();

    if (n <= size) {
        cout << "degree of f"<<n<<"is smaller than given size" << size << endl;
        return {};
    }

    vector<Fr> q(n - size);

    
    for (int i = n - 1; i >= size; i--) {
        Fr coeff = f[i];

        q[i - size] = coeff;

        f[i] -= coeff;       
        f[i - size] += coeff;   
    }

    for (int i = 0; i < size; i++) {
        if (!f[i].isZero()) {
            cout << "f is not divisible by X^size - 1" << endl;
            return {};
        }
    }

    return q;
}
//Prover
ZeroTestProof Prover_init(const vector<Fr>& f, int size, const Fr& r) {

    ZeroTestProof proof;

    PK pk = Setup((int)f.size() - 1);
    proof.pk = pk;
    proof.Cf = Commit(f, pk);

    // q(x) = f(x) / (x^l - 1)
    vector<Fr> q = vanishPoly(f, size);
    if (q.empty()) {
        cout << "ZeroTest prover failed." << endl;
        proof.valid = false;
        return proof; 
    }

    proof.Cq = Commit(q, pk);

    proof.Wf = CreateWitness(f, r, pk);
    proof.Wq = CreateWitness(q, r, pk);

    proof.valid = true;
    return proof;
}

//Verfier
bool ZeroTestVerifier(const ZeroTestProof& proof, int size, const Fr& r) {
    if (!proof.valid) {
        return false;
    }
    bool okF = VerifyEval(proof.pk, proof.Cf, proof.Wf);
    bool okQ = VerifyEval(proof.pk, proof.Cq, proof.Wq);

    if (!okF || !okQ) {
        return false;
    }
    if (proof.Wf.i != r || proof.Wq.i != r) {
        return false;
    }

    Fr rPowL;
    Fr fsize = size;
    Fr::pow(rPowL, r, fsize);

    Fr zh = rPowL - 1;

    Fr rhs = proof.Wq.value * zh;

    // Check f(r) == q(r)(r^l - 1)
    return proof.Wf.value == rhs;
}
void testZeroTest(const vector<Fr>& f, int size) {
    
    Fr r = random_fr();
    ZeroTestProof proof = Prover_init(f,size,r);
    bool ok = ZeroTestVerifier(proof, size, r);
    if (proof.valid) {
        cout << "ZeroTest example:" << endl;
        cout << "size = " << size << endl;
        cout << "r = " << r << endl;
        cout << "f(r) = " << proof.Wf.value << endl;
        cout << "q(r) = " << proof.Wq.value << endl;
        cout << "ZeroTest verification: "
            << (ok ? "passed" : "failed") << endl;
    }
}
int main() {
    init();
    // f(x) = -1 + x^4
    vector<Fr> f1 = {-1,0,0,0,1};
    vector<Fr> f2 = {-10,-9,-8,0,10,9,8};

    testZeroTest(f1,4);
    testZeroTest(f2,4);

    return 0;
}


#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <mcl/bn.hpp>

using namespace std;
using namespace mcl::bn;
using boost::multiprecision::cpp_int;

// BN254/BN_SNARK1 scalar field modulus
const string BN_SNARK1_FR_MOD =
    "21888242871839275222246405745257275088548364400416034343698204186575808495617";

// primitive root 
const string generator = "5";

cpp_int getMod() {
    return cpp_int(BN_SNARK1_FR_MOD);
}

Fr inverseFr(const Fr& x) {
    Fr temp;
    Fr::inv(temp, x);
    return temp;
}
mcl::Vint change(const cpp_int& x) {
    mcl::Vint v;
    v.setStr(x.convert_to<string>());
    return v;
}

bool rootOfUnity(Fr& omega, int n) {
    cpp_int r = getMod();

    if ((r - 1) % n != 0) {
        cout<<"n does not divide MOD - 1"<<endl;
        return false;
    }

    cpp_int exponent = (r - 1) / n;
    mcl::Vint exponentVint = change(exponent);

    Fr one = 1;

    for (int candidate = 2; candidate < 1000; candidate++) {
        Fr a;
        a = candidate;

        Fr tempomega;
        Fr::pow(tempomega, a, exponentVint);

        // Check omega^n == 1
        Fr temp1;
        Fr temp2;
        Fr::pow(temp1,tempomega, n);
        Fr::pow(temp2,tempomega, n/2);
        if (!(temp1 == one)) {
            continue;
        }

        if (n > 1 && (temp2 == one)) {
            continue;
        }

        //cout << "Candidate generator value = " << candidate << endl;
        omega = tempomega;
        return true;
    }
    
    cout<<"No valid root of unity found"<<endl;
    return false;
}

bool powerOfTwo(int n) {
    return n > 0 && ((n & (n - 1)) == 0);
}

void bitreverse(vector<Fr>& a) {
    int n = a.size();

    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        while (j & bit) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if (i < j) {
            swap(a[i], a[j]);
        }
    }
}

bool ntt(vector<Fr>& a, bool inv) {
    int n = a.size();
    if (!powerOfTwo(n)) {
        cout<<"NTT size must be power of 2"<<endl;
        return false;
    }

    cpp_int mod = getMod();
    if ((mod - 1) % n != 0) {
        cout<<"NTT size does not divide Fr modulus - 1"<<endl;
        return false;
    }
    bitreverse(a);

    for (int len = 2; len <= n; len <<= 1) {
        Fr wlen;
        if (!rootOfUnity(wlen, len)) {
            return false;
        }
        if (inv) {
            wlen = inverseFr(wlen);
        }
        for (int i = 0; i < n; i += len) {
            Fr w = 1;
            for (int j = 0;j < len / 2; j++) {
                Fr u= a[i + j];
                Fr v =a[i +j + len / 2] * w;
                a[i +j] = u + v;
                a[i+j +len / 2] = u - v;

                w *= wlen;
            }
        }
    }

    if (inv) {
        Fr nFr = n;
        Fr nInv = inverseFr(nFr);
        for (Fr& x : a) {
            x *= nInv;
        }
    }
    return true;
}

vector<Fr> randPoly(int n, int maxCoeff = 20) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dist(0, maxCoeff);

    vector<Fr> poly(n);
    for (int i = 0; i < n; i++) {
        poly[i] = dist(gen);
    }

    return poly;
}

bool equalPoly(const vector<Fr>& A, const vector<Fr>& B) {
    if (A.size() != B.size()) return false;
    for (int i = 0; i < A.size(); i++) {
        if (A[i] != B[i]) return false;
    }

    return true;
}
void printPoly(const string& label, const vector<Fr>& poly) {
    cout << label;
    for (const Fr& x : poly) {
        cout << x << " ";
    }
    cout << endl;
}
void testInter(const vector<Fr>& poly) {
    cout << "Testing NTT and INTT" << endl;
    cout << endl;

    printPoly("Original coefficients: ", poly);
    vector<Fr> values = poly;
    if (!ntt(values, false)) {
        cout << "NTT failed" << endl;
        return;
    }
    printPoly("NTT values: ", values);
    cout << endl;
    vector<Fr> res = values;
    if (!ntt(res, true)) {
        cout << "INTT failed" << endl;
        return;
    }
    printPoly("Recovered coefficients by INTT: ", res);

    cout << "Interpolation successful? ";
    cout << (equalPoly(poly, res) ? "YES" : "NO") << endl;
    cout << endl;
}

vector<Fr> polyMul(const vector<Fr>& A, const vector<Fr>& B) {
    int resultSize = A.size() + B.size() - 1;
    int n = 1;
    while (n < resultSize) {
        n <<= 1;
    }

    vector<Fr> FA(n);
    vector<Fr> FB(n);

    for (int i = 0; i < n; i++) {
        FA[i].clear();
        FB[i].clear();
    }
    for (int i = 0; i < A.size(); i++) {
        FA[i] = A[i];
    }
    for (int i = 0;i < B.size(); i++) {
        FB[i] = B[i];
    }
    if (!ntt(FA,false)) {
        cout << "NTT on A failed" << endl;
        return {};
    }
    if (!ntt(FB, false)) {
        cout << "NTT on B failed" << endl;
        return {};
    }
    vector<Fr> C(n);

    for (int i = 0; i < n; i++) {
        C[i] =FA[i]*FB[i];
    }

    if (!ntt(C, true)) {
        cout << "INTT on C failed" << endl;
        return {};
    }
    C.resize(resultSize);
    return C;
}

vector<Fr> normalPolyMul(const vector<Fr>& A, const vector<Fr>& B) {
    vector<Fr> C(A.size() + B.size() - 1);

    for (Fr& x : C) {
        x.clear();
    }
    for (int i = 0; i <A.size(); i++) {
        for (int j =0; j < B.size(); j++) {
            C[i + j] += A[i]*B[j];
        }
    }
    return C;
}

void testPolyMul(int a, int b) {
    cout << "Testing Polynomial Multiplication" << endl;
    cout << endl;

    vector<Fr> A = randPoly(a);
    vector<Fr> B = randPoly(b);

    vector<Fr> C_ntt = polyMul(A, B);
    vector<Fr> C_normal = normalPolyMul(A, B);

    printPoly("Poly A: ", A);
    printPoly("Poly B: ", B);
    printPoly("NTT multiplication result: ", C_ntt);
    printPoly("Normal multiplication result: ", C_normal);

    cout << "Multiplication correct? ";
    cout << (equalPoly(C_ntt, C_normal) ? "YES" : "NO") << endl;
    cout << endl;
}

int main() {
    initPairing(BN_SNARK1);

    int n = 16;

    vector<Fr> poly = randPoly(n);

    testInter(poly);
    int a = 5, b = 10;

    testPolyMul(a,b);

    return 0;
}
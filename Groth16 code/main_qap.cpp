#include "R1CS.hpp"
#include "QAP.hpp"

#include <iostream>
#include <vector>
#include <mcl/bn.hpp>
#include <random>

using namespace std;
using namespace mcl::bn;

int main() {
    initPairing(BN_SNARK1);

    // vector<Fr> coeff = {5,1,0,1};
    vector<Fr> coeff = rand_poly(5);
    print_poly(coeff);

    Fr x = 3;
    R1CS r1cs = poly_R1CS(coeff);
    vector<Fr> assign =witness_vector(coeff,x);
    cout << "Assignment:\n";
    for (int i = 0; i < assign.size(); i++) {
        cout << "z[" << i << "] = "<< assign[i] << endl;
    }

    bool valid =verify_R1CS(r1cs,assign);
    cout <<"\nPoly output y = " <<assign[1]<< endl;
    cout <<"Poly R1CS: "<<(valid ? "passed" :"failed")<< endl;

    QAP qap = convert(r1cs);
    print_checks(qap,assign);

    vector<Fr> h;
    bool valid_qap =verify_QAP(qap,assign,h);
    cout << "\nQAP verify: "<< (valid_qap ? "passed" : "failed")<< endl;
    cout << "\nh(X) coefficients:\n";
    for (int i=0;i< h.size();i++) {
        cout << "h[" << i << "] = "<< h[i]<< endl;
    }

    vector<Fr> wA = assign;
    wA[1] += 1;
    vector<Fr> wH;
    bool w_valid =verify_QAP(qap,wA,wH);
    cout << "\nWrong output test: "<< (w_valid? "incorrectly accepted": "correctly rejected")<< endl;
        
    return 0;
}
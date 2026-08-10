#include "R1CS.hpp"
#include "QAP.hpp"
#include "Groth16.hpp"

#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <mcl/bn.hpp>

using namespace std;
using namespace mcl::bn;
using namespace std::chrono;

int main() {
    initPairing(BN_SNARK1);
    // vector<Fr> coeff = {5,1,0,1};
    vector<Fr> coeff = rand_poly(9);
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

    
    Keys keys = groth_Setup(qap);

    cout << "\nGroth16 trusted setup completed."<< endl;
    cout << "Number of variables: "<< keys.pKey.num_v<< endl;
    cout << "Number of public inputs: "<< keys.pKey.PI<< endl;
    cout << "u-query G1 size: "<< keys.pKey.uQ_G1.size()<< endl;
    cout << "v-query G1 size: "<< keys.pKey.vQ_G1.size()<< endl;
    cout << "v-query G2 size: "<< keys.pKey.vQ_G2.size()<< endl;
    cout << "Public query size: "<< keys.vKey.public_G1.size()<< endl;
    cout << "Private query size: "<< keys.pKey.private_G1.size()<< endl;
    cout << "h-query size: "<< keys.pKey.hQ_G1.size()<< endl;

    GrothProof proof =groth_Prove(keys.pKey,qap,assign);
    cout << "\nGroth16 proof generated successfully."<< endl;

    vector<Fr> PI = {assign[1]};
    bool v =groth_Verify(keys.vKey,PI,proof);

    cout << "\nGroth16 proof verification: "<< (v ? "passed" : "failed")<< endl;

    vector<Fr> w = {assign[1] + Fr(1)};
    bool res2 =groth_Verify(keys.vKey,w,proof);
    cout << "Wrong public input test: "<< (res2? "incorrectly accepted": "correctly rejected")<< endl;
    

    return 0;
}
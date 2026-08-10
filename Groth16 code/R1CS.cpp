#include "R1CS.hpp"
#include <iostream>
#include <vector>
#include <random>

using namespace std;
using namespace mcl::bn;

vector<Fr> norm(vector<Fr> coeff) {
    if (coeff.empty()) {
        coeff.push_back(Fr(0));
    }
    while (coeff.size() > 1 && coeff.back().isZero()) {
        coeff.pop_back();
    }
    return coeff;
}
R1CS poly_R1CS(vector<Fr> coeff) {
    coeff = norm(coeff);

    int deg = static_cast<int>(coeff.size())- 1;

    R1CS r1cs;

    if(deg == 0){
        r1cs.v = 3;
    } else{
        r1cs.v = deg+2;
    }

    r1cs.PI = 1;


    if (deg == 0) {
        vector<Fr> rowA(r1cs.v,Fr(0));
        vector<Fr> rowB(r1cs.v,Fr(0));
        vector<Fr> rowC(r1cs.v,Fr(0));
        rowA[1] = 1;               
        rowA[0] -= coeff[0]; 
        rowB[0] = 1;          

        r1cs.A.push_back(rowA);
        r1cs.B.push_back(rowB);
        r1cs.C.push_back(rowC);
        return r1cs;
    }

    for (int i = deg-1;i >= 0;i--) {
        vector<Fr> rowA(r1cs.v,Fr(0));
        vector<Fr> rowB(r1cs.v, Fr(0));
        vector<Fr> rowC(r1cs.v, Fr(0));

        if (i + 1==deg) {
            rowA[0]=coeff[deg];
        } else {
            rowA[i+3] = 1;
        }

        rowB[2] = 1;

        if(i == 0){
            rowC[1] = 1;
        } else{
            rowC[i+2] = 1;
        }
        rowC[0] -= coeff[i];

        r1cs.A.push_back(rowA);
        r1cs.B.push_back(rowB);
        r1cs.C.push_back(rowC);
    }

    return r1cs;
}
vector<Fr> witness_vector(vector<Fr> coeff,const Fr& x) {
    coeff = norm(coeff);
    int deg = static_cast<int>(coeff.size()) - 1;

    int var;
    if(deg == 0){
        var = 3;
    } else{
        var = deg+2;
    }

    vector<Fr> assign(var, Fr(0));

    assign[0] = 1;
    assign[2] = x;
    if (deg == 0) {
        assign[1] = coeff[0];
        return assign;
    }

    vector<Fr> w(deg +1, Fr(0));

    w[deg] = coeff[deg];
    for (int i =deg-1;i>= 0;i--) {
        w[i] = w[i+1]*x+coeff[i];
    }

    assign[1] = w[0];
    for (int i = 1;i<deg;i++) {
        assign[i+2] = w[i];
        
    }

    return assign;
}
Fr multiply(const vector<Fr>& row, const vector<Fr>& assign) {
    Fr res = 0;

    for(int i=0;i<row.size();i++) {
        res += row[i]*assign[i];
    }
    return res;
}
bool verify_R1CS(const R1CS& r1cs,const vector<Fr>& assign) {
    if(assign.size()!=r1cs.v) {
        return false;
    }
    for(int i=0;i<r1cs.A.size();i++) {
        Fr l = multiply(r1cs.A[i], assign);
        Fr r = multiply(r1cs.B[i], assign);
        Fr out = multiply(r1cs.C[i], assign);

        cout << "\nConstraint " << i << ":" << endl;
        cout << "Left   = " << l << endl;
        cout << "Right  = " << r << endl;
        cout << "Output = " << out << endl;

        if (l*r !=out) {
            return false;
        }
    }
    return true;
}
Fr random_fr() {
    Fr x;
    x.setByCSPRNG();
    return x;
}
vector<Fr> rand_poly(int deg) {

    random_device rd;
    mt19937 gen(rd());

    uniform_int_distribution<int> dist(1, 30);

    vector<Fr> coeff(deg + 1);

    for (int i = 0; i <= deg; i++) {
        coeff[i] = dist(gen);
        //coeff[i] = random_fr();
    }

    return coeff;
}
void print_poly(const vector<Fr>& coeff) {
    cout << "Polynomial coefficients:" << endl;
    for (int i =0;i<coeff.size(); i++) {
        cout<<"c[" << i << "] = "<<coeff[i]<<endl;
    }
}

// int main() {
//     initPairing(BN_SNARK1);

//     // //f(x) = 5 + x + x^3
//     // vector<Fr> coeff = {5,1,0,1};
//     vector<Fr> coeff = rand_poly(5);
//     print_poly(coeff);

//     Fr x =3;
//     R1CS r1cs = poly_R1CS(coeff);
//     vector<Fr> a= witness_vector(coeff,x);
//     cout << "Assignment:"<< endl;

//     for (int i = 0; i < a.size(); i++) {
//         cout << "z["<< i<< "] = "
//              << a[i]<< endl;
//     }
//     bool valid =verify_R1CS(r1cs,a);
//     cout <<"\nPoly output y = " <<a[1]<< endl;
//     cout <<"Poly R1CS: "<<(valid ? "passed" :"failed")<< endl;

//     return 0;
// }
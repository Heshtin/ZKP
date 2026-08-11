#include "QAP.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>

using namespace std;
using namespace mcl::bn;

vector<Fr> norm_poly(vector<Fr> coeff) {
    if (coeff.empty()) {
        coeff.push_back(Fr(0));
    }
    while (coeff.size() > 1 && coeff.back().isZero()) {
        coeff.pop_back();
    }
    return coeff;
}

bool poly_zeroCheck(const vector<Fr>&poly) {
    for(const Fr& coeff:poly) {
        if(!coeff.isZero()){
            return false;
        }
    }
    return true;
}
Fr eval_poly(const vector<Fr>& poly,const Fr& x) {
    Fr res = 0;
    int n = static_cast<int>(poly.size()) - 1;
    for (int i= n;i>= 0;i--) {
        res*= x;
        res+= poly[i];
    }
    return res;
}

vector<Fr> add_poly(const vector<Fr>& f1,const vector<Fr>& f2) {
    int s = max(f1.size(), f2.size());
    vector<Fr> res(s, Fr(0));

    for (int i = 0; i<f1.size();i++) {
        res[i] += f1[i];
    }
    for (int i =0; i< f2.size();i++) {
        res[i] +=f2[i];
    }
    res = norm_poly(res);
    return res;
}

vector<Fr> sub_poly(const vector<Fr>& f1,const vector<Fr>& f2) {
    int s = max(f1.size(), f2.size());
    vector<Fr> res(s, Fr(0));

    for (int i = 0; i<f1.size();i++) {
        res[i] += f1[i];
    }
    for (int i =0; i<f2.size();i++) {
        res[i] -=f2[i];
    }
    return norm_poly(res);
}

vector<Fr> scale_poly(const vector<Fr>& poly,const Fr& s) {
    vector<Fr> res = poly;

    for (Fr& coeff : res) {
        coeff *= s;
    }
    return norm_poly(res);
}

vector<Fr> multi_poly(const vector<Fr>& f1,const vector<Fr>& f2) {
    if (f1.empty()||f2.empty()) {
        return {Fr(0)};
    }
    vector<Fr> res(f1.size()+f2.size()- 1,Fr(0));

    for (int i = 0; i < f1.size(); i++) {
        for (int j = 0; j < f2.size(); j++) {
            res[i+j]+= f1[i]*f2[j];
        }
    }
    return norm_poly(res);
}
vector<Fr> interpolate(const vector<Fr>& x,const vector<Fr>& y) {
    vector<Fr> res = {Fr(0)};
    for (int i = 0; i < x.size(); i++) {
        vector<Fr> basis = {Fr(1)};
        Fr d = 1;
        for (int k = 0; k<x.size();k++) {
            if (i==k) {
                continue;
            }
            vector<Fr> factor = {-x[k],Fr(1)};

            basis = multi_poly(basis,factor);

            d *= x[i] - x[k];
        }

        Fr inv_d;
        Fr::inv(inv_d, d);

        Fr mul =y[i]*inv_d;
        basis = scale_poly(basis,mul);

        res = add_poly(res,basis);
    }
    return norm_poly(res);
}
QAP convert(const R1CS& r1cs) {
    int num_c =static_cast<int>(r1cs.A.size());
    int num_v =r1cs.v;
    QAP qap;
    qap.num_c = num_c;
    qap.num_v = num_v;
    qap.PI =r1cs.PI;

    for (int i= 0; i<num_c; i++) {
        qap.pt.push_back(Fr(i+1));
    }
    qap.u.resize(num_v);
    qap.v.resize(num_v);
    qap.w.resize(num_v);

    for (int i=0;i<num_v;i++) {
        vector<Fr> aCol;
        vector<Fr> bCol;
        vector<Fr> cCol;
        for (int j= 0;j<num_c; j++) {
            aCol.push_back(r1cs.A[j][i]);
            bCol.push_back(r1cs.B[j][i]);
            cCol.push_back(r1cs.C[j][i]);
        }
        qap.u[i] = interpolate(qap.pt,aCol);
        qap.v[i] = interpolate(qap.pt,bCol);
        qap.w[i] = interpolate(qap.pt,cCol);
    }
    qap.t = {Fr(1)};

    for (const Fr& pt : qap.pt) {
        vector<Fr> factor = {-pt,Fr(1)};

        qap.t = multi_poly(qap.t,factor);
    }
    return qap;
}

vector<Fr> combine_poly(const vector<vector<Fr>>& poly,const vector<Fr>& assign) {
    vector<Fr> res = {Fr(0)};
    for (int i =0;i< assign.size();i++) {
        vector<Fr> term = scale_poly(poly[i],assign[i]);
        res = add_poly(res,term);
    }
    return res;
}
pair<vector<Fr>, vector<Fr>> divide_poly(vector<Fr> num,vector<Fr> deno) {
    num = norm_poly(num);
    deno= norm_poly(deno);
    if (poly_zeroCheck(deno)) {
        throw runtime_error("Division by zero");
    }

    if (poly_zeroCheck(num) ||num.size() < deno.size()) {
        return {{Fr(0)},num};
    }
    vector<Fr> q(num.size() - deno.size() + 1,Fr(0));

    Fr inv_coeff;

    Fr::inv(inv_coeff,deno.back());

    while (!poly_zeroCheck(num) && num.size() >= deno.size()) {
        int diff =num.size()-deno.size();
        Fr coeff =num.back()*inv_coeff;

        q[diff] += coeff;
        for (int i = 0; i<deno.size(); i++) {
            num[diff+ i] -=coeff* deno[i];
        }
        num = norm_poly(num);
    }

    q = norm_poly(q);
    num = norm_poly(num);
    return {q,num};
}

bool verify_QAP(const QAP& qap,const vector<Fr>& assign,vector<Fr>& h) {
    if (assign.size() !=qap.num_v) {
        return false;
    }
    vector<Fr> U =combine_poly(qap.u,assign);
    vector<Fr> V =combine_poly(qap.v,assign);
    vector<Fr> W =combine_poly(qap.w,assign);

    vector<Fr> nume =sub_poly(multi_poly(U, V),W);

    auto [q, r] =divide_poly(nume,qap.t);
    h = q;

    return poly_zeroCheck(r);
}

void print_checks(const QAP& qap,const vector<Fr>& assign) {
    vector<Fr> U =combine_poly(qap.u,assign);
    vector<Fr> V =combine_poly(qap.v,assign);
    vector<Fr> W =combine_poly(qap.w,assign);

    cout << "\nQAP constraint checks:\n";

    for (int i = 0;i <qap.pt.size(); i++) {
        const Fr& pt = qap.pt[i];

        Fr u_val = eval_poly(U, pt);
        Fr v_val= eval_poly(V, pt);
        Fr w_val = eval_poly(W, pt);

        cout << "\nConstraint " << i << endl;
        cout << "Point X = " << pt << endl;
        cout << " U(X) = "<< u_val
             << " V(X) = " << v_val
             << " W(X) = "<< w_val
             << endl;

        cout << "U(X)V(X) = "<< u_val * v_val<< endl;

        cout << (u_val * v_val == w_val? "passed": "failed") << endl;
    }
}
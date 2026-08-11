#include "Groth16.hpp"
#include "R1CS.hpp"
#include <iostream>
#include <stdexcept>

using namespace std;
using namespace mcl::bn;

Fr rand_fr() {
    Fr val;
    do {val.setByCSPRNG();
    } while (val.isZero());
    return val;
}
G1 multi_G1(const G1& pt,const Fr& s) {
    G1 res;
    G1::mul(res, pt, s);
    return res;
}
G2 multi_G2(const G2& pt,const Fr& s) {
    G2 res;
    G2::mul(res, pt, s);
    return res;
}

void addToG1(G1& acc,const G1& pt) {
    G1 res;
    G1::add(res, acc, pt);
    acc = res;
}

void addToG2(G2& acc,const G2& pt) {
    G2 res;
    G2::add(res, acc, pt);
    acc = res;
}

Keys groth_Setup(const QAP& qap) {

    Keys keys;

    G1 g1;
    G2 g2;
    mapToG1(g1, 1);
    mapToG2(g2, 1);

    Fr tau1;
    Fr tau2;

    do {
        tau1 = rand_fr();
        tau2 = eval_poly(qap.t, tau1);
    } while (tau2.isZero());

    Fr a = rand_fr();
    Fr b  = rand_fr();
    Fr g = rand_fr();
    Fr d = rand_fr();

    Fr g_inv;
    Fr d_inv;

    Fr::inv(g_inv, g);
    Fr::inv(d_inv, d);

    PKey& pk = keys.pKey;
    VKey& vk =keys.vKey;

    pk.a_G1 = multi_G1(g1, a);
    pk.b_G1 = multi_G1(g1, b);
    pk.b_G2 = multi_G2(g2, b);
    pk.d_G1 = multi_G1(g1, d);
    pk.d_G2 = multi_G2(g2, d);


    vk.a_G1 = pk.a_G1;
    vk.b_G2 = pk.b_G2;
    vk.g_G2 = multi_G2(g2, g);
    vk.d_G2 = pk.d_G2;

    pk.num_v =qap.num_v;
    pk.PI =qap.PI;
    vk.PI =qap.PI;
    for (int i = 0;i < qap.num_v;i++) {
        Fr ui_tau =eval_poly(qap.u[i], tau1);
        Fr vi_tau =eval_poly(qap.v[i], tau1);
        Fr wi_tau =eval_poly(qap.w[i], tau1);

        pk.uQ_G1.push_back(multi_G1(g1, ui_tau));

        pk.vQ_G1.push_back(multi_G1(g1, vi_tau));
        pk.vQ_G2.push_back(multi_G2(g2, vi_tau));

        Fr comb =b* ui_tau+ a *vi_tau+ wi_tau;
        if (i <= qap.PI) {
            Fr temp =comb *g_inv;
            vk.public_G1.push_back(multi_G1(g1,temp));
        } else {
            Fr temp = comb * d_inv;
            pk.private_G1.push_back(multi_G1(g1,temp));
        }
    }

    Fr tau_pow = 1;

    int max_deg =qap.num_c - 2;
    for (int i = 0;i <= max_deg;i++) {
        Fr s =tau_pow *tau2 *d_inv;
        pk.hQ_G1.push_back(multi_G1(g1, s));
        tau_pow *= tau1;
    }
    return keys;
}


GrothProof groth_Prove(const PKey& pk,const QAP& qap,const vector<Fr>& assign) {
    int num_v =static_cast<int>(assign.size());

    vector<Fr> h;

    bool valid_qap =verify_QAP(qap,assign,h);

    Fr r = rand_fr();
    Fr s = rand_fr();
    GrothProof proof;
    proof.A = pk.a_G1;
    proof.B = pk.b_G2;
    G1 b_in_G1 = pk.b_G1;

    for (int i= 0;i<num_v; i++) {
        G1 temp1 =multi_G1(pk.uQ_G1[i],assign[i]);
        addToG1(proof.A, temp1);

        G2 temp2 =multi_G2(pk.vQ_G2[i],assign[i]);
        addToG2(proof.B, temp2);

        G1 temp3 =multi_G1(pk.vQ_G1[i],assign[i]);
        addToG1(b_in_G1, temp3);
    }

    G1 r_G1 =multi_G1(pk.d_G1,r);
    addToG1(proof.A, r_G1);

    G2 s2_G2 =multi_G2(pk.d_G2,s);
    addToG2(proof.B, s2_G2);

    G1 s1_G1 =multi_G1(pk.d_G1,s);
    addToG1(b_in_G1, s1_G1);

    proof.C.clear();

    int f_PI =pk.PI + 1;

    for (int i = f_PI;i < num_v;i++) {
        int qi =i - f_PI;
        G1 temp =multi_G1(pk.private_G1[qi],assign[i]);
        addToG1(proof.C, temp);
    }

    for (int i = 0; i < h.size(); i++) {
        G1 temp =multi_G1(pk.hQ_G1[i],h[i]);
        addToG1(proof.C, temp);
    }

    G1 sA =multi_G1(proof.A,s);
    addToG1(proof.C, sA);

    G1 rB1 =multi_G1(b_in_G1,r);
    addToG1(proof.C, rB1);

    G1 rs =multi_G1(pk.d_G1,r * s);

    G1 finalC;
    G1::sub(finalC, proof.C, rs);
    proof.C = finalC;

    return proof;
}
bool groth_Verify(const VKey& vk,const vector<Fr>& PI,const GrothProof& proof) {

    if (PI.size() !=static_cast<size_t>(vk.PI)) {
        cout << "Incorrect number of public inputs."<< endl;
        return false;
    }

    if (vk.public_G1.size() !=PI.size() + 1) {
        cout << "Incorrect public-query size."<< endl;
        return false;
    }

    G1 public_acc =vk.public_G1[0];
    for (int i = 0; i <PI.size(); i++) {
        G1 temp = multi_G1(vk.public_G1[i+1],PI[i]);
        addToG1(public_acc,temp);
    }
    Fp12 lPair;

    pairing(lPair,proof.A,proof.B);

    Fp12 ab_pair;
    Fp12 public_pair;
    Fp12 c_pair;
    pairing(ab_pair,vk.a_G1,vk.b_G2);
    pairing(public_pair,public_acc,vk.g_G2);
    pairing(c_pair, proof.C, vk.d_G2);


    Fp12 rPair;
    Fp12::mul(rPair,ab_pair,public_pair);
    Fp12::mul(rPair,rPair,c_pair);

    return lPair == rPair;
}
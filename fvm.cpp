#include "geometry.hpp"
#include <iostream>
#include <math.h>
#include <string>
#include <list>
#include <tuple>

using namespace std;

// global index initialization for mesh 
const int nNode = 2000;
const int nEdge = 1000;
const int nCell = 800;

Node globalNodeList[nNode];
Edge globalEdgeList[nEdge];
Cell globalCellList[nCell];
list<int> NodeCellConnect[nNode];
list<int> CellEdgeConnect[nCell];
list<int> CellCellConnect[nCell];

float g = 9.81;
float pi = 3.1415926;
float delta_t = 0.1;


float sgn(float x){
    if (x>=0){
        return 1;
    }
    else{
        return 0;
    }
}


tuple<float, float, float> rotate(float u1, float u2, float u3, float theta){
    float u1_rot, u2_rot, u3_rot;
    u1_rot = u1;
    u2_rot = cosf(theta)*u2 + sinf(theta)*u3;
    u3_rot = -sinf(theta)*u2 + cosf(theta)*u3;
    return make_tuple(u1_rot, u2_rot, u3_rot);
}


// rotate a vector according to the rotation matrix
tuple<float, float, float> rotate(float u1, float u2, float u3, float theta){
    float u1_rot, u2_rot, u3_rot;
    u1_rot = u1;
    u2_rot = cosf(theta)*u2 + sinf(theta)*u3;
    u3_rot = -sinf(theta)*u2 + cosf(theta)*u3;
    return make_tuple(u1_rot, u2_rot, u3_rot);
}


// rotate a vector back according to the inverse of the rotation matrix
tuple<float, float, float> rotate_back(float u1, float u2, float u3, float theta){
    float u1_rot, u2_rot, u3_rot;
    u1_rot = u1;
    u2_rot = cosf(theta)*u2 - sinf(theta)*u3;
    u3_rot = sinf(theta)*u3 + cosf(theta)*u3;
    return make_tuple(u1_rot, u2_rot, u3_rot);
}

// Van Albada limiter for variable reconstruction
float VanAlbabaLimiter(float a, float b, float epsilon){
	if (a*b>0){
		return (powf(a, 2)*b+epsilon*b+powf(b, 2)*a+epsilon*a)/(powf(a, 2)+powf(b, 2)+2*epsilon);
	}
    else{
        return 0;
    }
}

// Roe average for U
float roeAverageU(float uL, float uR, float hL, float hR){
    float ur;
    if ((hL < 0) || (hR < 0)){
        cout << "Negative water depth at Roe average for u!" << endl;
    }
    if (hL*hL+hR*hR > 0){
        ur = (sqrtf(hL)*uL+sqrtf(hR)*uR)/(sqrtf(hL)+sqrtf(hR));
    }
    else{
        ur = 0;
    }
    return ur;
}

// Roe average for h
float roeAverageH(float hL, float hR){
    float hr;
    if ((hL < 0) || (hR < 0)){
        cout << "Negative water depth at Roe average for h!" << endl;
    }
    else{
        hr = sqrtf(hL*hR);
    }
    return hr;
}

// Roe average for c
float roeAverageC(float hL, float hR){
    float cr;
    if ((hL < 0) || (hR < 0)){
        cout << "Negative water depth at Roe average for c!" << endl;
    }
    else{
        cr = sqrtf(0.5*g*(hL+hR));
    }
    return cr;
}


// estimating nodal primitive varibles by Inverse Distance Weighting
void Node::u_estimation(){
    int l = NodeCellConnect[_Idx].size();
    float wsum = 0;
    float d = 0;
    // initilize conservative variables to 0
    _h = 0;
    _hu = 0;
    _hv = 0;
    // estimate by Inverse Distance Weighting
    list<int>::iterator i;
    for (i = NodeCellConnect[_Idx].begin(); i != NodeCellConnect[_Idx].end(); i++){
        d = sqrtf(powf(_px-globalCellList[*i].getValcpx(), 2)+powf(_py-globalCellList[*i].getValcpy(), 2));
        _h = _h + globalCellList[*i].getValh()/d;
        _hu = _hu + globalCellList[*i].getValhu()/d;
        _hv = _hv + globalCellList[*i].getValhv()/d;
        wsum += 1/d;
    }
    _h = _h/wsum;
    _hu = _hu/wsum;
    _hv = _hv/wsum;
}


void Edge::reconstruct(float epsilon_VAl, float epsilon_wd, float tol, float h0){
    if (boundaryFlag==0){
        reconstruct_MUSCL(epsilon_VAl, epsilon_wd);
    }
    else{
        reconstruct_FO();
        // for super-critical flow, u > sqrt(gh)
        if (cosf(theta_s)*_huL+sinf(theta_s)*_hvL >= _hL*sqrtf(g*_hL)){
            reconstruct_OpenBC_super();
        }
        // for sub-critical flow
        else if (boundaryFlag==2){
            reconstruct_OpenBC_h();
        }
        else if (boundaryFlag==3){
            reconstruct_OpenBC_u();
        }
        else if (boundaryFlag==4){
            reconstruct_OpenBC_hu(tol, h0);
        }
        // for solid boundary condition, we don't need the information at the right side
    }
}


void Edge::fluxCalc(){
    float S_mass, SL_mass, SR_mass;
    float d;

    d = min(AL, AR)/length;
    SL_mass = globalCellList[_owner].S;
    SR_mass = globalCellList[_neighbour].S;
    S_mass = 0.5*(SL_mass+SR_mass)*d;
    if (boundaryFlag==1){
        fluxCalc_SolidBC();
    }
    else{
        fluxCalc_hllcs(S_mass);
    }
}


void Edge::reconstruct_FO(){
    float fh, fhu, fhv, th, thu, thv;
    // get the variable value at the nodes of the edge
    fh = globalNodeList[_fromNode].getValh();
    fhu = globalNodeList[_fromNode].getValhu();
    fhv = globalNodeList[_fromNode].getValhv();
    th = globalNodeList[_toNode].getValh();
    thu = globalNodeList[_toNode].getValhu();
    thv = globalNodeList[_toNode].getValhv();
    // set variables value at the left side of the edge
    _hL = 0.5*(fh+th);
    _huL = 0.5*(fhu+thu);
    _hvL = 0.5*(fhv+thv);
}


// evaluate the left and right value of primitive variables at the midpoint
void Edge::reconstruct_MUSCL(float epsilon_VAl, float epsilon_wd){
    float mh, mhu, mhv;
    float fh, fhu, fhv, th, thu, thv;
    float oh, ohu, ohv, oz, nh, nhu, nhv, nz;
    float gradhLup, gradhuLup, gradhvLup;
    float gradhLdown, gradhuLdown, gradhvLdown;
    float gradhRup, gradhuRup, gradhvRup;
    float gradhRdown, gradhuRdown, gradhvRdown;
    // get the variable value at the nodes of the edge
    fh = globalNodeList[_fromNode].getValh();
    fhu = globalNodeList[_fromNode].getValhu();
    fhv = globalNodeList[_fromNode].getValhv();
    th = globalNodeList[_toNode].getValh();
    thu = globalNodeList[_toNode].getValhu();
    thv = globalNodeList[_toNode].getValhv();
    // get the variable value of the owner cell and the neighbour cell
    oh = globalCellList[_owner].getValh();
    ohu = globalCellList[_owner].getValhu();
    ohv = globalCellList[_owner].getValhv();
    nh = globalCellList[_neighbour].getValh();
    nhu = globalCellList[_neighbour].getValhu();
    nhv = globalCellList[_neighbour].getValhv();
    // variable value at the middle point of edge
    mh = 0.5*(fh+th);
    mhu = 0.5*(fhu+thu);
    mhv = 0.5*(fhv+thv);
    // Left side
    gradhLup = (oh - mh)/rcmL;
    gradhuLup = (ohu - mhu)/rcmL;
    gradhvLup = (ohv - mhv)/rcmL;
    gradhLdown = (globalNodeList[_pl].getValh() - oh)/rncL;
    gradhuLdown = (globalNodeList[_pl].getValhu() - ohu)/rncL;
    gradhvLdown = (globalNodeList[_pl].getValhv() - ohv)/rncL;
    // Right side
    gradhRup = (mh - nh)/rcmR;
    gradhuLup = (mhu - nhu)/rcmR;
    gradhvLup = (mhv - nhv)/rcmR;
    gradhLdown = (nh - globalNodeList[_pr].getValh())/rncR;
    gradhuLdown = (nhu - globalNodeList[_pr].getValhu())/rncR;
    gradhvLdown = (nhv - globalNodeList[_pr].getValhv())/rncR;
    // slope limiting
    _hL = oh + rcmL*VanAlbabaLimiter(gradhLdown, gradhLup, epsilon_VAl);
    _huL = ohu + rcmL*VanAlbabaLimiter(gradhuLdown, gradhuLup, epsilon_VAl);
    _hvL = ohv + rcmL*VanAlbabaLimiter(gradhvLdown, gradhvLup, epsilon_VAl);
    _hR = nh - rcmR*VanAlbabaLimiter(gradhRdown, gradhRup, epsilon_VAl);
    _huR = nhu - rcmR*VanAlbabaLimiter(gradhuRdown, gradhuRup, epsilon_VAl);
    _hvR = nhv - rcmR*VanAlbabaLimiter(gradhvRdown, gradhvRup, epsilon_VAl);
    // switch the second order MUSCL scheme to the first order scheme for wet-dry fronts
    if ((oh<=epsilon_wd) || (_hL<=min(fabsf(0.5*(zL-zR)), float(0.25)*oh))){
        _hL = oh;
        _huL = ohu;
        _hvL = ohv;
    }
    if ((nh<=epsilon_wd) || (_hR<=min(fabsf(0.5*(zL-zR)), float(0.25)*nh))){
        _hR = nh;
        _huR = nhu;
        _hvR = nhv;
    }
    // preserve non-negative water depth
    if (_hL<0){
        cout << "Reconstructed water depth is negative at the left side of edge " << _eIdx << endl;
        _hL = 0;
        _huL = 0;
        _hvL = 0;
    }
    if (_hR<0){
        cout << "Reconstructed water depth is negative at the right side of edge " << _eIdx << endl;
        _hR = 0;
        _huR = 0;
        _hvR = 0;
    }
}

// Boundary treatment
void Edge::reconstruct_OpenBC_h(){
    float uL_rot, vL_rot;
    float uR_rot, vR_rot;
    float uR_back, vR_back; 
    float hL_rot, huL_rot, hvL_rot;
    // get rotated u, v, i.e. normal velocity, tangential velocity
    tie(hL_rot, huL_rot, hvL_rot) = rotate(_hL, _huL, _hvL, theta_s);
    if (hL_rot>0){
        uL_rot = huL_rot/hL_rot;
        vL_rot = hvL_rot/hL_rot;
    }
    else {
        uL_rot = 0;
        vL_rot = 0;
    }

    uR_rot = uL_rot + 2*sqrtf(g)*(sqrtf(_hL)-sqrtf(_hL));
    vR_rot = vL_rot;
    uR_back = uR_rot*cosf(theta_s)-vR_rot*sinf(theta_s);
    vR_back = uR_rot*sinf(theta_s)+vR_rot*cosf(theta_s);
    _hR = _hBC;
    _huR = _hR*uR_back;
    _hvR = _hR*vR_back;
}

void Edge::reconstruct_OpenBC_u(){
    float uL_rot, vL_rot;
    float uR_rot, vR_rot;
    float uR_back, vR_back; 
    float hL_rot, huL_rot, hvL_rot;
    
    tie(hL_rot, huL_rot, hvL_rot) = rotate(_hL, _huL, _hvL, theta_s);
    if (hL_rot>0){
        uL_rot = huL_rot/hL_rot;
        vL_rot = hvL_rot/hL_rot;
    }
    else {
        uL_rot = 0;
        vL_rot = 0;
    }

    vR_rot = vL_rot;
    uR_back = uR_rot*cosf(theta_s)-vR_rot*sinf(theta_s);
    vR_back = uR_rot*sinf(theta_s)+vR_rot*cosf(theta_s);
    _hR = powf(uL_rot+2*sqrtf(g*hL_rot)-_uBC, 2)/(4*g);
    _huR = _hR*uR_back;
    _hvR = _hR*vR_back;
}

void Edge::reconstruct_OpenBC_hu(float tol, float h0){
    float uL_rot, vL_rot;
    float uR_rot, vR_rot;
    float uR_back, vR_back; 
    float hL_rot, huL_rot, hvL_rot;
    float C, res, h_sqrt;
    
    tie(hL_rot, huL_rot, hvL_rot) = rotate(_hL, _huL, _hvL, theta_s);
    if (hL_rot>0){
        uL_rot = huL_rot/hL_rot;
        vL_rot = hvL_rot/hL_rot;
    }
    else {
        uL_rot = 0;
        vL_rot = 0;
    }

    // solve hR_rot, uR_rot by:
    // 1) hR_rot * uR_rot = huBC
    // 2) uR_rot + 2*sqrt(g*hR_rot) = uL_rot + 2*sqrt(g*hL_rot)
    // Using Newton-Raphson method
    C = uL_rot + 2*sqrtf(g*hL_rot);
    h_sqrt = sqrtf(h0);
    res = 2*sqrtf(g)*h_sqrt*h_sqrt*h_sqrt - C*h_sqrt*h_sqrt + _huBC;
    while (res > tol){
        h_sqrt -= res/(6*sqrtf(g)-2*C*h_sqrt);
        res = 2*sqrtf(g)*h_sqrt*h_sqrt*h_sqrt - C*h_sqrt*h_sqrt + _huBC;
    }
    uR_rot = _huBC/(h_sqrt*h_sqrt);
    vR_rot = vL_rot;
    
    uR_back = uR_rot*cosf(theta_s)-vR_rot*sinf(theta_s);
    vR_back = uR_rot*sinf(theta_s)+vR_rot*cosf(theta_s);
    _hR = h_sqrt*h_sqrt;
    _huR = _hR*uR_back;
    _hvR = _hR*vR_back;
}

void Edge::reconstruct_OpenBC_super(){
    _hR = _hL;
    _huR = _huL;
    _hvR = _huR;
}

// HLLC solver implementation
// Make sure both sides of one water depth is positive
void Edge::fluxCalc_hllcs(float S_mass){
    float Sl, Sr, Ss_p, Ss_n;
    float h_roe, c_roe, u_roe;
    float hL_rot, huL_rot, hvL_rot;
    float uL_rot, vL_rot;
    float hR_rot, huR_rot, hvR_rot;
    float uR_rot, vR_rot;
    float d;
    float S, S_z, S_tau;
    float dL, dR, delta_d, delta_z, delta_zi, T, T_max;
    float H1;
    float h_hllcs, hu_hllcs, hv_hllcs;
    float phi1L_rot, phi2L_rot, phi3L_rot;
    float phi1R_rot, phi2R_rot, phi3R_rot;
    // rotate the primitive variable vector; flux should be pointed to neighbour
    tie(hL_rot, huL_rot, hvL_rot) = rotate(_hL, _huL, _hvL, theta_s);
    tie(hR_rot, huR_rot, hvR_rot) = rotate(_hR, _huR, _hvR, theta_s);
    if (hL_rot>0){
        uL_rot = huL_rot/hL_rot;
        vL_rot = hvL_rot/hL_rot;
    }
    else {
        uL_rot = 0;
        vL_rot = 0;
    }
    if (hR_rot>0){
        uR_rot = huR_rot/hR_rot;
        vR_rot = hvR_rot/hR_rot;
    }
    else {
        uR_rot = 0;
        vR_rot = 0;
    }
    // left wave and right wave speed estimation
    c_roe = roeAverageC(hL_rot, hR_rot);
    if (hL_rot*hL_rot+hR_rot*hR_rot > 0){
        u_roe = roeAverageU(uL_rot, uR_rot, hL_rot, hR_rot);
    }
    else {
        u_roe = 0;
    }
    Sl = u_roe - c_roe;
    Sr = u_roe + c_roe;
    // source terms integral
    // 1) bed slope term
    dL = zL + hL_rot;
    dR = zR + hR_rot;
    delta_d = dR - dL;
    delta_z = zR - zL;
    // -- determine delta_zi
    if ((delta_z>=0) and (dL<zR)){
        delta_zi = hL_rot;
    }
    else if ((delta_z<0) and (dR<zL)){
        delta_zi = hR_rot;
    }
    else{
        delta_zi = delta_z;
    }
    // -- determine T
    if (delta_z>=0){
        T = g*(hL_rot-fabsf(delta_zi)/2)*delta_zi;
    }
    else{
        T = g*(hR_rot-fabsf(delta_zi)/2)*delta_zi;
    }
    h_roe = roeAverageH(hL_rot, hR_rot);
    if (fabsf(g*h_roe*delta_z)>=fabsf(T)){
        T_max = g*h_roe*delta_z;
    }
    else{
        T_max = T;
    }
    if ((delta_d*delta_z>=0) and (u_roe*delta_z)>0){
        S_z = T_max;
    }
    else{
        S_z = T;
    }
    // 2) fricition term
    d = min(AL, AR)/length;
    if (h_roe>0){
        // fricition term is estimated by (Delta x)*Cf*u*|u|,
        // where an empirical estimation for Cf is Cf = g*n^2/h^{1/3}
        S_tau = min(d*g*n*n/powf(h_roe, 1/3)*u_roe*fabsf(u_roe), sgn(u_roe)*fabsf(h_roe*fabsf(u_roe)*min(fabsf(uL_rot), fabsf(uR_rot))));
    }
    else{
        S_tau = 0;
    }
    S = S_z + S_tau;
    // check if Sl >=0 or Sr <= 0 
    if (Sl>=0){
        phi1L_rot = huL_rot;
        phi1R_rot = phi1L_rot + S_mass;
        phi2L_rot = huL_rot*uL_rot + 1/2*g*hL_rot*hL_rot;
        phi2R_rot = phi2L_rot + S;
        phi3L_rot = hvL_rot*uL_rot;
        phi3R_rot = phi3L_rot;
    }
    else if (Sr<=0){
        phi1R_rot = huR_rot;
        phi1L_rot = phi1R_rot - S_mass;
        phi2R_rot = huR_rot*uR_rot + 1/2*g*hR_rot*hR_rot;
        phi2L_rot = phi2R_rot - S;
        phi3R_rot = hvR_rot*uR_rot;
        phi3L_rot = phi3R_rot;
    }
    // if not, estimate middle wave speed for HLLCS solver
    else{
        H1 = 2*u_roe/(u_roe*u_roe-c_roe*c_roe)*S_mass-1/(u_roe*u_roe-c_roe*c_roe)*S;
        Ss_p = (Sl*hR_rot*(uR_rot-Sr)-Sr*hL_rot*(uL_rot-Sl)+Sr*Sl*H1-Sr*S_mass)/(hR_rot*(uR_rot-Sr)-hL_rot*(uL_rot-Sl)+Sl*H1-S_mass);
        Ss_p = (Ss_p+fabsf(Ss_p))/2;
        Ss_n = (Sl*hR_rot*(uR_rot-Sr)-Sr*hL_rot*(uL_rot-Sl)+Sr*Sl*H1-Sl*S_mass)/(hR_rot*(uR_rot-Sr)-hL_rot*(uL_rot-Sl)+Sr*H1-S_mass);
        Ss_n = (Ss_n-fabsf(Ss_p))/2;
        if (Ss_p*Ss_n == 0){
            // positive middle wave speed
            if (Ss_p > 0){
                h_hllcs = (hL_rot*(uL_rot-Sl)-Sl*H1+S_mass)/(Ss_p-Sl);
                hu_hllcs = h_hllcs*Ss_p;
                hv_hllcs = h_hllcs*vL_rot;
                phi1L_rot = huL_rot + Sl*(h_hllcs-hL_rot);
                phi1R_rot = phi1L_rot + S_mass;
                phi2L_rot = huL_rot*uL_rot + 1/2*g*hL_rot*hL_rot + Sl*(hu_hllcs-huL_rot);
                phi2R_rot = phi2L_rot + S;
                phi3L_rot = hvL_rot*uL_rot + Sl*(hv_hllcs-hvL_rot);
                phi3R_rot = phi3L_rot;
            }
            // negative middle wave speed
            else{
                h_hllcs = (hR_rot*(uR_rot-Sr)+Sr*H1-S_mass)/(Ss_n-Sr);
                hu_hllcs = h_hllcs*Ss_n;
                hv_hllcs = h_hllcs*vR_rot;
                phi1R_rot = huR_rot + Sr*(h_hllcs-hR_rot);
                phi1L_rot = phi1R_rot - S_mass;
                phi2R_rot = huR_rot*uR_rot + 1/2*g*hR_rot*hR_rot + Sr*(hu_hllcs-huR_rot);
                phi2L_rot = phi2R_rot - S;
                phi3R_rot = hvR_rot*uR_rot + Sr*(hv_hllcs-hvR_rot);
                phi3L_rot = phi3R_rot;
            }
        }
        // assuming that the advection terms are dominant with respect to the source terms
        else if (sgn(Ss_p) == sgn(u_roe)){
            h_hllcs = (hL_rot*(uL_rot-Sl)-Sl*H1+S_mass)/(Ss_p-Sl);
            hu_hllcs = h_hllcs*Ss_p;
            hv_hllcs = h_hllcs*vL_rot;
            phi1L_rot = huL_rot + Sl*(h_hllcs-hL_rot);
            phi1R_rot = phi1L_rot + S_mass;
            phi2L_rot = huL_rot*uL_rot + 1/2*g*hL_rot*hL_rot + Sl*(hu_hllcs-huL_rot);
            phi2R_rot = phi2L_rot + S;
            phi3L_rot = hvL_rot*uL_rot + Sl*(hv_hllcs-hvL_rot);
            phi3R_rot = phi3L_rot;
        }
        else{
            h_hllcs = (hR_rot*(uR_rot-Sr)+Sr*H1-S_mass)/(Ss_n-Sr);
            hu_hllcs = h_hllcs*Ss_n;
            hv_hllcs = h_hllcs*vR_rot;
            phi1R_rot = huR_rot + Sr*(h_hllcs-hR_rot);
            phi1L_rot = phi1R_rot - S_mass;
            phi2R_rot = huR_rot*uR_rot + 1/2*g*hR_rot*hR_rot + Sr*(hu_hllcs-huR_rot);
            phi2L_rot = phi2R_rot - S;
            phi3R_rot = hvR_rot*uR_rot + Sr*(hv_hllcs-hvR_rot);
            phi3L_rot = phi3R_rot;
        }
    }
    tie(_phi1L, _phi2L, _phi3L) = rotate_back(phi1L_rot, phi2L_rot, phi3L_rot, theta_s);
    tie(_phi1R, _phi2R, _phi3R) = rotate_back(phi1R_rot, phi2R_rot, phi3R_rot, theta_s);
}


void Edge::fluxCalc_SolidBC(){
    _phi1R = _phi2R = _phi3R = 0;
    _phi1L = 0;
    _phi2L = 0.5*g*_hL*_hL*cosf(theta_s);
    _phi3L = 0.5*g*_hL*_hL*sinf(theta_s);
}


void Cell::step(){
    list<int>::iterator j;
    float phi1, phi2, phi3, l;
    // get the node of cell
    for (j = CellEdgeConnect[_cIdx].begin(); j != CellEdgeConnect[_cIdx].end(); j++){
        if (globalEdgeList[*j].getValOCell()==_cIdx){
            l = globalEdgeList[*j].length;
            phi1 = globalEdgeList[*j].getValPhi1L();
            phi2 = globalEdgeList[*j].getValPhi2L();
            phi3 = globalEdgeList[*j].getValPhi3L();
            _h = _h - phi1*l/area*delta_t;
            _hu = _hu - phi2*l/area*delta_t;
            _hv = _hv - phi3*l/area*delta_t;
        }
        else {
            phi1 = globalEdgeList[*j].getValPhi1R();
            phi2 = globalEdgeList[*j].getValPhi2R();
            phi3 = globalEdgeList[*j].getValPhi3R();
            _h = _h + phi1*l/area*delta_t;
            _hu = _hu + phi2*l/area*delta_t;
            _hv = _hv + phi3*l/area*delta_t;
        }
    }
}
#include "geometry.hpp"
#include <math.h>
#include <list>
#include <set>
#include <iostream>
#include <algorithm>

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

// constructor for Node
Node::Node(float x, float y, float ele, int i){
    _Idx = i;
    _px = x;
    _py = y;
    _b = ele;
}

Node::Node(const Node &rhs){
    _Idx = rhs.getValId();
    _px = rhs.getValpx();
    _py = rhs.getValpy();
    _b = rhs.getValb();
}

// constructor for Edge
Edge::Edge(int eIdx, int fId, int tId, int flag){
    _eIdx = eIdx;
    _fromNode = fId;
    _toNode = tId;
    _owner = -1;  
	_neighbour = -1; 
    _hL = _hR = _huL = _huR = _hvL = _hvR = 0;
    _phi1L = _phi2L = _phi3L = _phi1R = _phi2R = _phi3R = 0;
    _hBC = _uBC = _huBC;
    zL = zR = 0;
    AL = AR = 0;
    n = 0;
    rcmL = rncL = rcmR = rncR = 0;
    // boundaryFlag = 0: inner edge
    // boundaryFlag = 1: solid wall boundary
    // boundaryFlag = 2: open boundary where h is given
    // boundaryFlag = 3: open boundary where uN is given
    // boundaryFlag = 4: open boundary where h*uN is given
    boundaryFlag = flag; 
    length = sqrtf(powf(globalNodeList[fId].getValpx()-globalNodeList[tId].getValpx(), 2)+powf(globalNodeList[fId].getValpy()-globalNodeList[tId].getValpy(), 2));
    if (globalNodeList[fId].getValpx() == globalNodeList[tId].getValpx()){
        theta = pi/2;
    }
    else{
        theta = atanf((globalNodeList[fId].getValpy()-globalNodeList[tId].getValpy()/(globalNodeList[fId].getValpx()-globalNodeList[tId].getValpx())));
    }
    theta_s = pi/2 + theta;
}

// constructor for Cell
Cell::Cell(int cIdx, int e1, int e2, int e3, float n_r, float h, float hu, float hv){
    int N1, N2, N3;
    _cIdx = cIdx;
    _h = h;
    _hu = hu;
    _hv = hv;
    n = n_r;
    N1 = globalEdgeList[e1].getValFrom();
    N2 = globalEdgeList[e1].getValTo();
    if ((globalEdgeList[e2].getValFrom() == N1) || (globalEdgeList[e2].getValFrom() == N2)){
        N3 = globalEdgeList[e2].getValFrom();
    }
    else{
        N3 = globalEdgeList[e2].getValTo();
    }
    // center elevation and coordinates
    bc = (globalNodeList[N1].getValb()+globalNodeList[N2].getValb()+globalNodeList[N3].getValb())/3;
    _cpx = (globalNodeList[N1].getValpx()+globalNodeList[N2].getValpx()+globalNodeList[N3].getValpx())/3;
    _cpy = (globalNodeList[N1].getValpy()+globalNodeList[N2].getValpy()+globalNodeList[N3].getValpy())/3;
    // calculate area
    area = ((globalNodeList[N2].getValpx()-globalNodeList[N1].getValpx())*(globalNodeList[N3].getValpy()-globalNodeList[N1].getValpy())-
            (globalNodeList[N3].getValpx()-globalNodeList[N1].getValpx())*(globalNodeList[N2].getValpy()-globalNodeList[N1].getValpy()))/2;
    area = fabsf(area);
    // calculate chi
    chi = 0;
    chi = min(area/globalEdgeList[e1].length, area/globalEdgeList[e2].length, area/globalEdgeList[e3].length);
    // set this cell as the OwnerCell for nodes
    NodeCellConnect[N1].push_back(cIdx);
    NodeCellConnect[N2].push_back(cIdx);
    NodeCellConnect[N3].push_back(cIdx);
}

Cell::Cell(const Cell &rhs){
    _cIdx = rhs._cIdx;
    _h = rhs.getValh();
    _hu = rhs.getValhu();
    _hv = rhs.getValhv();
}

// Given Cell-Edge connection, construct topological relationship for Edge
int topSearch(Node *gNlist, Edge *gElist, Cell *gClist, list<int> *CEconnect){
    int f, t, o, n, nId;
    float mx, my, rcm;
    float cmx, cmy, z, A, n;
    // determine the topological relationship of the edge
    for (int i=0; i<nCell; i++){
        Cell ci = gClist[i]; 
        int cid = ci.getValId();
        set<int> Ncset;
        set<int> Neset;  // nodes of edge
        set<int> Noset;  // opposite node
        list<int>::iterator j;
        // get the node of cell
        for (j = CEconnect[cid].begin(); j != CEconnect[cid].end(); j++){
            f = gElist[*j].getValFrom();
            t = gElist[*j].getValTo();
            if (Ncset.count(f) == 0){
                Ncset.insert(f);
            }
            if (Ncset.count(t) == 0){
                Ncset.insert(t);
            }
        }
        // set the owner cell, neightbour cell and the opposite points of the edge
        for (j=CEconnect[cid].begin(); j!=CEconnect[cid].end(); i++){
            f = gElist[*j].getValFrom();
            t = gElist[*j].getValTo();
            Neset.insert(f);
            Neset.insert(t);
            mx = (gNlist[f].getValpx()+gNlist[t].getValpx())/2;
            my = (gNlist[f].getValpy()+gNlist[t].getValpy())/2;
            cmx = gClist[cid].getValcpx()-mx;
            cmy = gClist[cid].getValcpy()-my;
            z = gClist[cid].bc;
            A = gClist[cid].area;
            n = gClist[cid].n;
            // if the owner cell of this edge has not been set
            // then set this cell as its owner cell and also set its left opposite point
            set_difference(Ncset.begin(), Ncset.end(), Neset.begin(), Neset.end(), inserter(Noset, Noset.end()));
            nId = *Noset.begin();
            rcm = sqrtf(powf(cmx, 2)+powf(cmy, 2));
            if (gElist[*j].getValOCell() < 0){
                gElist[*j].setValOCell(cid);
                gElist[*j].setValPl(nId);
                gElist[*j].zL = z;
                gElist[*j].AL = A;
                gElist[*j].n = n;
                gElist[*j].rcmL = rcm;
                gElist[*j].rncL = 2*rcm;
                if (cosf(gElist[*j].theta_s)*cmx+sinf(gElist[*j].theta_s)*cmy < 0){
                    gElist[*j].theta_s += pi;
                }
            }
            // else the owner cell of this edge has been set
            // then set this cell as its neighbour cell and also set its right opposite point
            else{
                gElist[*j].setValNCell(ci.getValId());
                gElist[*j].setValPr(nId);
                gElist[*j].zR = z;
                gElist[*j].AR = A;
                gElist[*j].n = (n+gElist[*j].n)/2;
                gElist[*j].rcmR = rcm;
                gElist[*j].rncR = 2*rcm;
            }
            Neset.clear();
            Noset.clear();
        }
    }
    // determine the neighbour cell of a cell
    for (int i=0; i<nEdge; i++){
        o = gElist[i].getValOCell();
        n = gElist[i].getValNCell();
        // For boundary egdes, n = -1
        if (n >= 0){
            CellCellConnect[n].push_back(o);
            CellCellConnect[o].push_back(n);
        }
        else{
            gElist[i].boundaryFlag = true;
        }
    }
}
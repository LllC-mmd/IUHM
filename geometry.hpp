#include <vector>

using namespace std;

class Node{
	private:
		int _Idx;
		float _b;
		float _px;
		float _py;
		// nodal primitive variable value, [h, hu, hv]
		float _h;
		float _hu;
		float _hv; 
	public:
		// Node constructor
		Node();
		Node(float x, float y, float ele, int i); 
		Node(const Node &rhs);  // copy constructor
		// get and set value
		int getValId() const {return _Idx;};
		float getValb() const {return _b;};
		float getValpx() const {return _px;};
		float getValpy() const {return _py;};
		float getValh() const {return _h;};
		float getValhu() const {return _hu;};
		float getValhv() const {return _hv;};
		void setValh(float h){_h=h;};
		void setValhu(float hu){_hu=hu;};
		void setValhv(float hv){_hv=hv;};
		// estimate nodal value for U by Inverse Distance Interpolation
		void u_estimation();
		// Destructor
		~Node();
};

class Edge{
	private:
		int _eIdx;
		int _fromNode;
		int _toNode;
		// index of the owner cell and the neighbor cell
		int _owner;  
		int _neighbour; 
		// index of the node across the edge, from left/owner to right/neighbor
		int _pl;
		int _pr;
		float _hL;
		float _huL;
		float _hvL;
		float _hR;
		float _huR;
		float _hvR;
		float _phi1L;
		float _phi1R;
		float _phi2L;
		float _phi2R;
		float _phi3L;
		float _phi3R;
		// boundary condition
		float _hBC;
		float _uBC;
		float _huBC;
	public:
		float theta;
		float theta_s;
		float length;
		float zL;
		float zR;
		float AL;
		float AR;
		float n;
		float rcmL;  // distance between the central point of the cell and the mid point of the edge
    	float rncL;  // distance between the central point of the cell and the left node across the edge
    	float rcmR;
    	float rncR;
		int boundaryFlag;
		// Edge constructor
		Edge();
		Edge(int eIdx, int fId, int tId, int flag);
		// get and set value
		int getValId(){return _eIdx;};
		int getValFrom(){return _fromNode;};
		int getValTo(){return _toNode;};
		int getValOCell(){return _owner;};
		int getValNCell(){return _neighbour;};
		float getValhL(){return _hL;};
		float getValhuL(){return _huL;};
		float getValhvL(){return _hvL;};
		float getValhR(){return _hR;};
		float getValhuR(){return _huR;};
		float getValhvR(){return _hvR;};
		float getValPhi1L(){return _phi1L;};
		float getValPhi1R(){return _phi1R;};
		float getValPhi2L(){return _phi2L;};
		float getValPhi2R(){return _phi2R;};
		float getValPhi3L(){return _phi3L;};
		float getValPhi3R(){return _phi3R;};
		void setValhL(float hL){_hL=hL;};
		void setValhuL(float huL){_huL=huL;};
		void setValhvL(float hvL){_hvL=hvL;};
		void setValhR(float hR){_hR=hR;};
		void setValhuL(float huR){_huR=huR;};
		void setValhvL(float hvR){_hvR=hvR;};
		void setValOCell(int p){_owner=p;};
		void setValNCell(int p){_neighbour=p;};
		void setValPl(int p){_pl=p;};
		void setValPr(int p){_pr=p;};
		// Finite Volume Method
		void reconstruct(float epsilon_VAl, float epsilon_wd, float tol, float h0);  // variable reconstruction at the mid point
		void reconstruct_MUSCL(float epsilon_VAl, float epsilon_wd); // reconstruction using MUSCL
		void reconstruct_FO(); // reconstruction using the first order scheme for the left side of boundary edges
		void reconstruct_OpenBC_h();  // reconstruction at the open boundary where h is given
		void reconstruct_OpenBC_u();  // reconstruction at the open boundary where u is given
		void reconstruct_OpenBC_hu(float tol, float h0);  // reconstruction at the open boundary where Q=hu is given
		void reconstruct_OpenBC_super(); // reconstruction at the open boundary for supercritical flow
		void fluxCalc();
		void fluxCalc_hllcs(float S_mass);  // determine the flux caused by source terms using HLLCS solver
		void fluxCalc_SolidBC();  // For solid boundary, flux can be directly computed
		// Destructor
		~Edge();
};

class Cell{
	private:
		int _cIdx;
		float _cpx;
		float _cpy;
		float _h;
		float _hu;
		float _hv;
	public:
		// physical property
		float S; // source term caused by effective rainfall, gully drainage
		int label;
		float area;
		float n;
		float bc;  // centroid bed elevation
		float zeta; // water surface level above reference surface, h = zeta - bc
		float d;  // infiltration, inlet flow to sewer or river, p - d = pe
		float slx;   // slope along the x direction
		float sly;   // slope along the y direction
		float chi;  // relevant distance used in the HLLCS solver
		// Cell constructor
		Cell();
		Cell(int cIdx, int e1, int e2, int e3, float n_r, float h, float hu, float hv);
		Cell(const Cell &rhs);
		// get and set value
		int getValId() const {return _cIdx;};
		float getValcpx() const {return _cpx;};
		float getValcpy() const {return _cpy;};
		float getValh() const {return _h;};
		float getValhu() const {return _hu;};
		float getValhv() const {return _hv;};
		void setValh(float h){_h=h;};
		void setValhu(float hu){_hu=hu;};
		void setValhv(float hv){_hv=hv;};
		void setValS(float S_mass){S=S_mass;};
		// Finite Volume Method
		void step();  // use fluxes of edges to update the cell-averaged variable
		// int inundation_check();  
		// Destructor
		~Cell();
};

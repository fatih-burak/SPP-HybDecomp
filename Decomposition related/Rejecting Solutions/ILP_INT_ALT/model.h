#ifndef BM_H
#define BM_H
using namespace std;
#include "gurobi_c++.h"
#include "helper_functions.h"
#include "slave_functions.h"

void ILP_INT_ALT(Instance& inst, Solution& sol);

class mycallbackT2 : public GRBCallback
{
public:
	Instance& inst;  // Store a reference to inst
	Solution& sol;
	const vector<vector<GRBVar>>& x;  // existing member variables
	const int H;
	const vector<vector<bool>>& NPs;

	vector<vector<int>> solution;
	// Constructor
	mycallbackT2(Instance& instS, Solution& solS, const vector<vector<GRBVar>>& xS, const int HS, const vector<vector<bool>>& NPsS)
		: inst(instS), sol(solS), x(xS), H(HS), NPs(NPsS) {}  // Initialize inst as a reference
protected:
	void callback();
};

#endif 


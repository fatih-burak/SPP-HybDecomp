#ifndef BM_H
#define BM_H



using namespace std;
#include "gurobi_c++.h"
#include "helper_functions.h"
#include "slave_functions.h"


void MP_ILP_DEL_LCBC(Instance& inst, Solution& sol);

class mycallbackT1 : public GRBCallback
{
public:
	Instance inst;
	Solution sol;
	vector<vector<GRBVar>> x;
	int H;

	mycallbackT1(const Instance& instS, Solution& sol, const vector<vector<GRBVar> >& x);

protected:
	void callback();
};



#endif 


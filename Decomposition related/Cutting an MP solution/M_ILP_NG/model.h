#ifndef MODEL_H
#define MODEL_H

using namespace std;
#include <iostream>


#include "gurobi_c++.h"
#include "helper_functions.h"
#include "slave_functions.h"


void ILP_NG(Instance& inst, Solution& sol);

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


#ifndef LAZYCALLBACKSPP_H_INCLUDED
#define LAZYCALLBACKSPP_H_INCLUDED
#include "model.h"
/****************************************************
 *  CLASS LazyCallbackSPP
 **********************/
 /*
 class LazyCallbackSPP : public IloCplex::LazyConstraintCallbackI 
{
public:
	Instance& inst;  // Store a reference to inst
	Solution& sol;
	const vector<vector<IloIntVar>>& x;  // existing member variables
	const int H;
	const vector<vector<bool>>& NPs;

	vector<vector<int>> solution;
    // Constructor
    LazyCallbackSPP(IloEnv env, Instance& instS, Solution& solS, const vector<vector<IloIntVar>>& xS, const int HS, const vector<vector<bool>>& NPsS);  // Initialize inst as a reference;

    IloCplex::CallbackI* duplicateCallback() const override {
        return new (getEnv()) LazyCallbackSPP(*this);
    }

    void main() override;
};

*/
/****************************************************
 *  CLASS LazyCallbackSPP
 **********************/
class IncumbentCallbackSPP : public IloCplex::IncumbentCallbackI
{
public:
	Instance& inst;  // Store a reference to inst
	Solution& sol;
	const vector<vector<IloIntVar>>& x;  // existing member variables
	const int H;
	const vector<vector<bool>>& NPs;

	vector<vector<int>> solution;
	// Constructor
	IncumbentCallbackSPP(IloEnv env, Instance& instS, Solution& solS, const vector<vector<IloIntVar>>& xS, const int HS, const vector<vector<bool>>& NPsS);  // Initialize inst as a reference;

	IloCplex::CallbackI* duplicateCallback() const override {
		return new (getEnv()) IncumbentCallbackSPP(*this);
	}

	void main() override;
};

#endif // LAZYCALLBACK2C_H_INCLUDED

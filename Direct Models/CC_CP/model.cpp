#include "model.h"
#include <random>

void CC_CP(Instance& inst, Solution& sol) {
	int start = getCPUTime();

	// from CSP to BPP
	vector<vector<int>> items2;
	int n2 = 0;
	for (int j = 0; j < inst.n; j++) {
		for (int k = 0; k < inst.items[j][2]; k++) {
			items2.push_back({ inst.items[j][0],inst.items[j][1],j });
			cout << "item=" << n2 << ": " << inst.items[j][0] << ", " << inst.items[j][1] << ", " << j << endl;
			n2++;
		}
	}
	inst.items2 = items2;

	// create a model
	IloEnv env;
	IloModel model(env);

	// define IloIntExprArray to store the end points of interval variables
	IloIntExprArray ends(env);
	// declaration of the variables for the model
	IloIntervalVarArray x(env, n2);
	IloIntervalVarArray y(env, n2);
	IloCumulFunctionExpr z(env);
	IloCumulFunctionExpr q(env);

	// initizalization of the variables for the model
	for (int j = 0; j < n2; j++) {
		x[j] = IloIntervalVar(env);
		x[j].setStartMin(0);
		x[j].setStartMax(inst.W - items2[j][0]);
		x[j].setSizeMin(items2[j][0]);
		x[j].setSizeMax(items2[j][0]);
		z += IloPulse(x[j], items2[j][1]);

		y[j] = IloIntervalVar(env);
		y[j].setStartMin(0);;
		ends.add(IloEndOf(y[j]));
		y[j].setSizeMin(items2[j][1]);
		y[j].setSizeMax(items2[j][1]);
		q += IloPulse(y[j], items2[j][0]);
	}

	// set the cumul constraint
	model.add(q <= inst.W);

	//Non-overlap constraints
	for (int j = 0; j < n2; j++) {
		for (int i = j + 1; i < n2; i++) {
			IloOr OR_const(env);
			OR_const.add(IloEndOf(x[i]) <= IloStartOf(x[j]));
			OR_const.add(IloEndOf(x[j]) <= IloStartOf(x[i]));
			OR_const.add(IloEndOf(y[i]) <= IloStartOf(y[j]));
			OR_const.add(IloEndOf(y[j]) <= IloStartOf(y[i]));
			model.add(OR_const);
		}
	}


	// set the objective: minimize the makespan
	IloIntVar obj(env, 0, IloInfinity); //you can set a LB on z
	model.add(obj >= IloMax(ends));
	model.add(obj >= z); // adding it makes the model faster
	model.add(obj >= sol.iLB);

	//minimize objective
	IloObjective objective = IloMinimize(env, obj);
	model.add(objective);

	// change some settings
	IloCP cp(model);
	cp.setParameter(IloCP::Workers, 1);
	cp.setParameter(IloCP::LogPeriod, 500000);
	cp.setParameter(IloCP::FailureDirectedSearch, IloCP::On);
	cp.setParameter(IloCP::CumulFunctionInferenceLevel, IloCP::Extended);
	cp.setParameter(IloCP::SearchType, IloCP::Restart);
	cp.setParameter(IloCP::TimeLimit, inst.timeLimit);
	cp.setParameter(IloCP::RandomSeed, 5);
	cp.propagate();
	//SOLVE
	cp.solve();

	int status = cp.getStatus();
	sol.UB = 9999 - sol.prepacked_height;
	sol.LB = cp.getObjBound();

	if (status == 2) { // optimal solution is found
		sol.UB = cp.getObjValue();
		cout << "Optimal found" << endl; sol.opt = 1;
		sol.LB = sol.UB;
		for (int j = 0; j < n2; j++) {
			cout << "Rectangle-" << j << " of size=" << items2[j][0] << "x" << items2[j][1] << " x coordinate = " << cp.getStart(x[j]) << ",y coordinate = " << cp.getStart(y[j]) << endl;
		}
	}
	else if (status == 0 || (status == 1 && (getCPUTime() - start + 2) > inst.timeLimit)) { //time limit reached
		cout << "Time limit is reached." << endl;
		if (status == 1) {
			cout << "However, some feasible solutions are found (which are not proven optimal)." << endl;
			sol.UB = cp.getObjValue();
		}
	}
	else if (status == 3) { //model is infeasible
		cout << "Model is infeasible." << endl;
		sol.opt = -1;
	}

	cp.end(); env.end();
	return;
}
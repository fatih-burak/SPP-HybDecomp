#include "model.h"
#include <random>

void DC_CP(Instance& inst, Solution& sol) {
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

	// local variable
	vector<bool> isWActive(inst.W, false);

	// compute normal patterms
	vector<vector<bool> > NPs(n2, vector<bool>(inst.W, false));
	for (int j = 0; j < n2; j++) {
		NPs[j][0] = true;
		for (int jp = 0; jp < n2; jp++) {
			if (jp != j) {
				for (int i = inst.W - items2[jp][0] - 1; i >= 0; i--) {
					if (NPs[j][i]) NPs[j][i + items2[jp][0]] = true;
				}
			}
		}
	}
	inst.NPs = NPs;

	// create a model
	IloEnv env;
	IloModel model(env);

	// define IloIntExprArray to store the end points of interval variables
	IloIntExprArray ends(env);
	// create & initialize y_modes[j][p] variables
	IloArray<IloIntervalVarArray> y_modes(env, n2);

	IloSolution startSol(env); // ---- WARM START ----
	for (int j = 0; j < n2; j++) {
		int p_warm = get<0>(sol.packing_locations[j]);// ---- WARM START ----
		int q_warm = get<1>(sol.packing_locations[j]);// ---- WARM START ----

		y_modes[j] = IloIntervalVarArray(env, inst.W - items2[j][0] + 1);
		for (int p = 0; p <= inst.W - items2[j][0]; p++) {
			if (NPs[j][p]) {
				isWActive[p] = true;
				y_modes[j][p] = IloIntervalVar(env);
				y_modes[j][p].setOptional();
				y_modes[j][p].setStartMin(0);
				y_modes[j][p].setSizeMin(items2[j][1]);
				y_modes[j][p].setSizeMax(items2[j][1]);
				ends.add(IloEndOf(y_modes[j][p]));
				if (p == p_warm) { // ---- WARM START ----
					startSol.setStart(y_modes[j][p], q_warm);// ---- WARM START ----
					startSol.setPresent(y_modes[j][p]);// ---- WARM START ----
				}// ---- WARM START ----
			}
		}

		

		IloExpr modeSelected(env);
		for (int p = 0; p <= inst.W - items2[j][0]; p++) {
			if (NPs[j][p]) modeSelected += IloPresenceOf(env, y_modes[j][p]);
		}
		model.add(modeSelected == 1);  // One mode must be selected
		modeSelected.end();
	}

	IloArray<IloIntervalVarArray> items_in_column(env, inst.W);
	for (int q = 0; q < inst.W; q++) {
		items_in_column[q] = IloIntervalVarArray(env);
	}

	for (int j = 0; j < n2; j++) {
		for (int q = 0; q <= inst.W - items2[j][0]; q++) {
			for (int p = q; p < q + items2[j][0]; p++) {
				if (NPs[j][q]) items_in_column[p].add(y_modes[j][q]);
			}
		}
	}

	// no-overlap constraints
	for (int q = 0; q < inst.W; q++) {
		if (isWActive[q]) model.add(IloNoOverlap(env, items_in_column[q]));
	}

	// set the objective: minimize the makespan
	IloIntVar z(env, 0, IloIntMax);
	model.add(z >= IloMax(ends));
	IloObjective objective = IloMinimize(env, z);
	model.add(objective);
	//bound z
	//model.add(z >= sol.iLB);

	// change some settings
	IloCP cp(model);
	cp.setStartingPoint(startSol);// ---- WARM START ----

	cp.setParameter(IloCP::Workers, 1);
	cp.setParameter(IloCP::LogPeriod, 500000);
	//cp.setParameter(IloCP::FailureDirectedSearch, IloCP::On);
	cp.setParameter(IloCP::CumulFunctionInferenceLevel, IloCP::Extended);
	cp.setParameter(IloCP::SearchType, IloCP::Restart);
	cp.setParameter(IloCP::TimeLimit, inst.timeLimit);
	cp.propagate();
	//SOLVE
	cp.solve();

	int status = cp.getStatus();
	sol.LB = cp.getObjBound();

	if (status == 2) { // optimal solution is found
		sol.UB = cp.getObjValue();
		cout << "Optimal found" << endl; sol.opt = 1;
		sol.LB = sol.UB;
		for (int j = 0; j < n2; j++) {
			for (int q = 0; q <= inst.W - items2[j][0]; q++) {
				if (NPs[j][q] && cp.isPresent(y_modes[j][q])) {
					cout << "Rectangle-" << j << " of size=" << items2[j][0] << "x" << items2[j][1] << " x coordinate = " << q << ",y coordinate = " << cp.getStart(y_modes[j][q]) << endl;
					break;
				}
			}
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
#include "model.h"


void D_ILP(Instance& inst, Solution& sol) {
	double startModel = getCPUTime();
	// from CSP to BPP
	vector<vector<int> > items2;
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
	vector<bool> isHActive(sol.UB, false);

	// compute normal patterns: width
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

	// compute normal patterns: height
	vector<vector<bool> > H_NPs(n2, vector<bool>(sol.UB, false));
	for (int j = 0; j < n2; j++) {
		H_NPs[j][0] = true;
		for (int jp = 0; jp < n2; jp++) {
			if (jp != j) {
				for (int i = sol.UB - items2[jp][1] - 1; i >= 0; i--) {
					if (H_NPs[j][i]) H_NPs[j][i + items2[jp][1]] = true;
				}
			}
		}
	}

	inst.NPs = NPs;
	inst.H_NPs = H_NPs;

	// create a model
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 1);
	GRBModel model = GRBModel(env);

	// declaration of the variables for the model
	GRBVar z = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
	GRBLinExpr obj = z;

	vector<vector<vector<GRBVar>>> x;
	// Resize the 3D vector
	x.resize(inst.W, vector<vector<GRBVar>>(sol.UB, vector<GRBVar>(n2)));

	// initizalization of the variables for the model
	for (int i = 0; i < n2; i++) {
		int p_warm = get<0>(sol.packing_locations[i]); 	// ---- WARM START ----
		int q_warm = get<1>(sol.packing_locations[i]); 	// ---- WARM START ----

		for (int p = 0; p <= inst.W - items2[i][0]; p++) {
			for (int q = 0; q <= sol.UB - items2[i][1]; q++) {
				if (NPs[i][p] && H_NPs[i][q]) {
					isWActive[p] = true;
					isHActive[q] = true;
					x[p][q][i] = model.addVar(0, 1, 0, GRB_BINARY);
					if (p == p_warm && q == q_warm) x[p][q][i].set(GRB_DoubleAttr_Start, 1.0);// ---- WARM START ----
					else x[p][q][i].set(GRB_DoubleAttr_Start, 0.0);// ---- WARM START ----
				}
			}
		}
	}
	// create linear expressions
	vector<GRBLinExpr> assigned(n2, 0);
	vector<vector<GRBLinExpr>> non_overlap;
	non_overlap.resize(inst.W, vector<GRBLinExpr>(sol.UB));
	vector<GRBLinExpr> height2(inst.W, 0); //extra constraint

	for (int i = 0; i < n2; i++) {
		for (int p = 0; p <= inst.W - items2[i][0]; p++) {
			for (int q = 0; q <= sol.UB - items2[i][1]; q++) {
				if (NPs[i][p] && H_NPs[i][q]) {
					assigned[i] += x[p][q][i];
					for (int pp = p; pp < p + items2[i][0]; pp++) {
						height2[pp] += items2[i][1] * x[p][q][i];
						for (int qq = q; qq < q + items2[i][1]; qq++) {
							non_overlap[pp][qq] += x[p][q][i];
						}
					}
				}
			}
		}
	}

	vector<vector<GRBLinExpr>> height;
	int start_h = sol.iLB;
	height.resize(n2, vector<GRBLinExpr>(sol.UB));
	for (int i = 0; i < n2; i++) {
		if (start_h <= sol.iLB - items2[i][1]) start_h = sol.iLB - items2[i][1];
		for (int q = sol.iLB - items2[i][1]; q <= sol.UB - items2[i][1]; q++) {
			for (int p = 0; p <= inst.W - items2[i][0]; p++) {
				if (NPs[i][p] && H_NPs[i][q]) {
					height[i][q] += (q + items2[i][1]) * x[p][q][i];
				}
			}
		}
	}

	//constraint set to bound z
	for (int i = 0; i < n2; i++) {
		for (int q = start_h; q <= sol.UB - items2[i][1]; q++) {
			model.addConstr(height[i][q] <= z);
		}
	}

	// create assignment constraints
	for (int i = 0; i < n2; i++) {
		model.addConstr(assigned[i] == 1);
	}

	// create height constraints 
	for (int p = 0; p < inst.W; p++) {
		if (isWActive[p]) {
			model.addConstr(height2[p] <= z);
			for (int q = 0; q < sol.UB; q++) {
				if (isHActive[q]) {
					model.addConstr(non_overlap[p][q] <= 1);
				}
			}
		}
	}
	model.update();
	// set the objective: minimize z
	model.setObjective(obj, GRB_MINIMIZE);

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_DoubleParam_TimeLimit, inst.timeLimit); //use the remaining time left if you increase the LB
	//model.getEnv().set(GRB_DoubleParam_BestObjStop, sol.iLB + 0.01); //activate for LB models

	// optimize	
	model.optimize();

	// store the results in a Solution object
	sol.Nvar = model.get(GRB_IntAttr_NumVars);
	sol.Nconstr = model.get(GRB_IntAttr_NumConstrs);
	sol.Ncoeff = model.get(GRB_IntAttr_NumNZs);

	sol.UB = 9999 - sol.prepacked_height;
	sol.LB = ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);
	if (model.get(GRB_IntAttr_SolCount) >= 1) {
		sol.UB = ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);
		if (sol.LB == sol.UB || sol.iLB == sol.UB) {
			//if (sol.LB == sol.UB) {
			sol.opt = 1;
			if (sol.LB != sol.UB) cout << "Optimal found but LB is not increased by the solver; instead BestObjStop is used." << endl;
		}
	}

	if (model.get(GRB_IntAttr_Status) == 3) { //Infeasible
		cout << "Infeasible..." << endl;
		sol.opt = -1;
	}
	if (model.get(GRB_IntAttr_Status) == 9) {
		cout << "Time limit is reached." << endl;
	}
	sol.timeT += getCPUTime() - startModel;
}


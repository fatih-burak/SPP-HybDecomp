#include "model.h"

void DC_ILP(Instance& inst, Solution& sol) {
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
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	// declaration of the variables for the model
	vector<vector<GRBVar> > x;
	GRBVar z = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
	GRBLinExpr obj = z;
	x.resize(inst.W, vector<GRBVar>(n2));

	// initizalization of the variables for the model
	for (int j = 0; j < n2; j++) {
		for (int i = 0; i <= inst.W - items2[j][0]; i++) {
			if (NPs[j][i]) {
				isWActive[i] = true;
				x[i][j] = model.addVar(0, 1, 0, GRB_BINARY);
			}
		}
	}

	vector<GRBVar> y(n2);
	for (int j = 0; j < n2; j++) {
		y[j] = model.addVar(0, sol.UB - items2[j][1], 0, GRB_CONTINUOUS);
	}
	
	// ---- WARM START ----
	for (int i = 0; i < n2; i++) {
		int p = get<0>(sol.packing_locations[i]);
		int q = get<1>(sol.packing_locations[i]);
		if (!NPs[i][p]) {
			cout << "Not in normal pattern!" << endl;
			continue;
		}
		x[p][i].set(GRB_DoubleAttr_Start, 1.0);
		y[i].set(GRB_DoubleAttr_Start, q);
		for (int pp = 0; pp <= inst.W - items2[i][0]; pp++) {
			if (NPs[i][pp] && pp!=p) {
				x[pp][i].set(GRB_DoubleAttr_Start, 0.0);
			}
		}
	}

	// declaration of the variables for the model
	vector<vector<GRBVar> > b;
	b.resize(n2, vector<GRBVar>(n2));
	for (int i = 0; i < n2; i++) {
		for (int j = i + 1; j < n2; j++) {
			b[i][j] = model.addVar(0, 1, 0, GRB_BINARY);
		}
	}

	// create linear expressions
	vector<GRBLinExpr> assigned(n2, 0);
	vector<GRBLinExpr> height(inst.W, 0);
	vector<GRBLinExpr> y_height(n2, 0);

	vector<vector<GRBLinExpr>> engages;
	engages.resize(inst.W, vector<GRBLinExpr>(n2, 0));

	for (int j = 0; j < n2; j++) {
		y_height[j] += y[j] + items2[j][1];
		for (int i = 0; i <= inst.W - items2[j][0]; i++) {
			if (NPs[j][i]) {
				assigned[j] += x[i][j];
				for (int ip = i; ip < i + items2[j][0]; ip++) {
					height[ip] += items2[j][1] * x[i][j];
					engages[ip][j] += x[i][j];
				}
			}
		}
	}

	// create assignment constraints
	for (int j = 0; j < n2; j++) {
		model.addConstr(assigned[j] == 1);
	}

	// create height constraints (pulse) 
	for (int i = 0; i < inst.W; i++) {
		if (isWActive[i])
			model.addConstr(height[i] <= z);
	}

	// create height constraints (pulse) 
	for (int j = 0; j < n2; j++) {
		model.addConstr(y_height[j] <= z);
	}

	int M = sol.UB;
	// Add non-overlap constraints
	for (int i = 0; i < n2; i++) {
		for (int j = i + 1; j < n2; j++) {
			for (int q = 0; q < inst.W; q++) {
				if (isWActive[q] = false) continue;

				// Constraint 1: y_i + h_i <= y_j + M(1 - b_ij) + M(2 - sum(x_ip) - sum(x_jq))
				model.addConstr(y[i] + items2[i][1] <= y[j] + M * (1 - b[i][j]) + M * (2 - engages[q][i] - engages[q][j]),
					"nonoverlap1_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(q));

				// Constraint 2: y_j + h_j <= y_i + M*b_ij + M(2 - sum(x_ip) - sum(x_jq))
				model.addConstr(y[j] + items2[j][1] <= y[i] + M * b[i][j] + M * (2 - engages[q][i] - engages[q][j]),
					"nonoverlap2_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(q));
			}
		}
	}
	// set the objective: minimize z
	model.setObjective(obj, GRB_MINIMIZE);

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_DoubleParam_TimeLimit, inst.timeLimit); //use the remaining time left if you increase the LB
	model.getEnv().set(GRB_DoubleParam_BestObjStop, sol.iLB + 0.01);
	// optimize	
	model.optimize();

	// store the results in a Solution object
	sol.Nvar = model.get(GRB_IntAttr_NumVars);
	sol.Nconstr = model.get(GRB_IntAttr_NumConstrs);
	sol.Ncoeff = model.get(GRB_IntAttr_NumNZs);

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
}

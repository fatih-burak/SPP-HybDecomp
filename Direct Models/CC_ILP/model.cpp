#include "model.h"

void CC_ILP(Instance& inst, Solution& sol) {
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
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	// declaration of the variables for the model
	GRBVar z = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
	GRBLinExpr obj = z;

	vector<GRBVar> x(n2), y(n2);
	for (int i = 0; i < n2; i++) {
		x[i] = model.addVar(0, inst.W - items2[i][0], 0, GRB_CONTINUOUS);
		y[i] = model.addVar(0, sol.UB - items2[i][1], 0, GRB_CONTINUOUS);
	}

	// ---- WARM START ----
	for (int i = 0; i < n2; i++) {
		int p = get<0>(sol.packing_locations[i]);
		int q = get<1>(sol.packing_locations[i]);
		x[i].set(GRB_DoubleAttr_Start, p);
		y[i].set(GRB_DoubleAttr_Start, q);
	}

	vector<vector<GRBVar>> h, v; //horizontal position, vertical position
	h.resize(n2, vector<GRBVar>(n2)); v.resize(n2, vector<GRBVar>(n2));
	for (int i = 0; i < n2; i++) {
		for (int j = 0; j < n2; j++) {
			if (i == j) continue;
			h[i][j] = model.addVar(0, 1, 0, GRB_BINARY);
			v[i][j] = model.addVar(0, 1, 0, GRB_BINARY);
		}
	}

	// create height bounding constraints
	for (int i = 0; i < n2; i++) model.addConstr(y[i] + items2[i][1] <= z);

	// create nonoverlap constraints
	for (int i = 0; i < n2; i++) {
		for (int j = i + 1; j < n2; j++) {
			// create OR selection for nonoverlap
			model.addConstr(h[i][j] + h[j][i] + v[i][j] + v[j][i] == 1);
			//horizontal nonoverlap
			model.addConstr(x[i] + items2[i][0] <= x[j] + inst.W * (1 - h[i][j]));
			model.addConstr(x[j] + items2[j][0] <= x[i] + inst.W * (1 - h[j][i]));
			//vertical nonoverlap
			model.addConstr(y[i] + items2[i][1] <= y[j] + sol.UB * (1 - v[j][i]));
			model.addConstr(y[j] + items2[j][1] <= y[i] + sol.UB * (1 - v[i][j]));
		}
	}

	model.update();
	// set the objective: minimize z
	model.setObjective(obj, GRB_MINIMIZE);

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_DoubleParam_TimeLimit, inst.timeLimit); 
	//model.getEnv().set(GRB_DoubleParam_BestObjStop, sol.iLB + 0.01); //activate for LB models
	// optimize	
	model.optimize();

	// store the results in a Solution object
	sol.Nvar = model.get(GRB_IntAttr_NumVars);
	sol.Nconstr = model.get(GRB_IntAttr_NumConstrs);
	sol.Ncoeff = model.get(GRB_IntAttr_NumNZs);

	sol.LB = ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);
	if (model.get(GRB_IntAttr_SolCount) >= 1) {
		sol.UB = ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);
		//if (sol.LB == sol.UB || sol.iLB == sol.UB) { //for LB models activate
		if (sol.LB == sol.UB) {
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

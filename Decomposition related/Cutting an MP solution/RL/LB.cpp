#include "LB.h"


void LB(const Instance& inst, Solution& sol) {
	// find the LB solution
	int LB_1 = LB_CBP_CONF(inst, sol);
	int LB_2 = LB_RE(inst, sol);
	int LB_3 = LB_PCC_CONF(inst, sol);
	cout << "LB_CBP_CONF=" << LB_1 << " LB RE=" << LB_2 << " LB_PCC_CONF=" << LB_3 << endl;
	sol.iLB += max({ LB_1, LB_2, LB_3 });
	sol.LB = sol.iLB;
	sol.vLB1 += sol.prepacked_height; sol.vLB2 += sol.prepacked_height; sol.vLB3 += sol.prepacked_height; //update the LBs by adding the prepacked height
}

int LB_CBP_CONF(const Instance& inst, Solution& sol) {
	double start = getCPUTime();
	// create node matrix
	vector<vector<int> > matToId(inst.items.size() + 1, vector<int>(inst.W + 1, -1));
	vector<vector<int> > idToMat;
	int nodeCount = 0;

	// create item arcs
	vector<bool> isA(inst.W + 1, false); isA[0] = true;
	vector<vector<int> > arcs;
	int tail, head;
	for (int j = 0; j < inst.items.size(); j++) {
		for (int i = inst.W; i >= 0; i--) {
			if (isA[i]) {
				if (matToId[j][i] == -1) {
					matToId[j][i] = nodeCount;
					idToMat.push_back({ j,i });
					nodeCount++;
				}
				tail = matToId[j][i];
				for (int k = 0; k <= inst.items[j][2]; k++) {
					if (i + inst.items[j][0] * k <= inst.W) {
						if (matToId[j + 1][i + inst.items[j][0] * k] == -1) {
							matToId[j + 1][i + inst.items[j][0] * k] = nodeCount;
							idToMat.push_back({ j + 1,i + inst.items[j][0] * k });
							nodeCount++;
						}
						head = matToId[j + 1][i + inst.items[j][0] * k];
						arcs.push_back({ tail,head,j,k });
						isA[i + inst.items[j][0] * k] = true;
					}
					else
						break;
				}
			}
		}
	}

	// Create a Gurobi environment

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	model.getEnv().set(GRB_IntParam_OutputFlag, 0); //silence the output

	// create a model           
	GRBLinExpr obj = 0;

	// declaration of the variables for the model
	vector<GRBVar> x(arcs.size());

	// initizalization of the variables for the model
	for (int k = 0; k < arcs.size(); k++)
		x[k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);

	model.update();

	// create linear expressions
	vector<GRBLinExpr> assigned(inst.items.size(), 0);
	vector<GRBLinExpr> cIn(idToMat.size(), 0);
	vector<GRBLinExpr> cOut(idToMat.size(), 0); if (cOut.size() == 0) cOut.resize(1, 0);
	for (int k = 0; k < arcs.size(); k++) {
		assigned[arcs[k][2]] += x[k] * arcs[k][3];
		cIn[arcs[k][1]] += x[k];
		cOut[arcs[k][0]] += x[k];
	}

	model.update();

	// create assignment constraints
	for (int j = 0; j < inst.items.size(); j++)
		model.addConstr(assigned[j] == inst.items[j][2] * inst.items[j][1]);

	// create flow conservation constraints 
	for (int i = 1; i < idToMat.size(); i++)
		model.addConstr(cIn[i] >= cOut[i]);

	// set the objective: minimize z
	model.setObjective(cOut[0], GRB_MINIMIZE);

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_IntParam_Method, 2);
	model.getEnv().set(GRB_IntParam_MIPFocus, 1);
	model.getEnv().set(GRB_DoubleParam_TimeLimit, 1);

	// find the optimal solution		
	model.optimize();

	try {
		sol.tLB += getCPUTime() - start;
		sol.tLB1 = getCPUTime() - start;
		sol.vLB1 = (int)ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);
		return (int)ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);

	}
	catch (GRBException e) {
		cout << "time limit in LB calculation" << endl;
		//std::cerr << "Gurobi Exception: " << e.getMessage() << std::endl;
		//std::cerr << "Error code: " << e.getErrorCode() << std::endl;
		sol.vLB2 = 0;
		return 0;
	}
	catch (...) {
		//std::cerr << "An unexpected error occurred." << std::endl;
		sol.vLB2 = 0;
		return 0;
	}

}


int LB_RE(const Instance& inst, Solution& sol) {
	double start = getCPUTime();
	vector<vector<int> > items2 = inst.items;

	// scale the instance
	for (int j = 0; j < inst.n; j++)
		items2[j][0] *= 2;

	// compute normal patterms
	vector<vector<bool> > NPs(inst.n, vector<bool>(inst.W, false));
	for (int j = 0; j < inst.n; j++) {
		NPs[j][0] = true;
		for (int jp = 0; jp < inst.n; jp++) {
			int lim = items2[jp][2];
			if (jp == j) lim--;
			for (int l = 0; l < lim; l++) {
				for (int i = inst.W - items2[jp][0] - 1; i >= 0; i--) {
					if (NPs[j][i]) NPs[j][i + items2[jp][0]] = true;
				}
			}
		}
	}

	// Create a Gurobi environment
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	model.getEnv().set(GRB_IntParam_OutputFlag, 0); //silence the output

	// declaration of the variables for the model
	vector<vector<GRBVar> > x;
	x.resize(inst.W, vector<GRBVar>(inst.n));
	//GRBVar z = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	GRBVar z = model.addVar(sol.vLB1, GRB_INFINITY, 0, GRB_CONTINUOUS); // start from the best known LB

	GRBLinExpr obj = z;

	// initizalization of the variables for the model
	for (int j = 0; j < inst.n; j++) {
		for (int i = 0; i < inst.W; i++) {
			if (NPs[j][i] && i <= inst.W - abs(inst.W - (i + items2[j][0]))) {
				x[i][j] = model.addVar(0, items2[j][2], 0, GRB_CONTINUOUS);
			}
		}
	}
	model.update();

	// create linear expressions
	vector<GRBLinExpr> assigned(inst.n, 0);
	vector<GRBLinExpr> height(inst.W, 0);
	for (int j = 0; j < inst.n; j++) {
		for (int i = 0; i < inst.W; i++) {
			if (NPs[j][i] && i <= inst.W - abs(inst.W - (i + items2[j][0]))) {
				assigned[j] += x[i][j];
				for (int ip = i; ip < min(inst.W, i + items2[j][0]); ip++) {
					height[ip] += items2[j][1] * x[i][j];
				}
				for (int ip = inst.W - 1; ip >= 2 * inst.W - (i + items2[j][0]); ip--) {
					height[ip] += items2[j][1] * x[i][j];
				}
			}
		}
	}
	model.update();

	// create assignment constraints
	for (int j = 0; j < inst.n; j++)
		model.addConstr(assigned[j] == items2[j][2]);

	// create height constraints 
	for (int i = 0; i < inst.W; i++) {
		model.addConstr(height[i] <= 2 * z);
	}

	// set the objective: minimize z
	model.setObjective(obj, GRB_MINIMIZE);

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_DoubleParam_TimeLimit, 1);

	// find the optimal solution		
	model.optimize();

	sol.tLB += getCPUTime() - start;
	sol.tLB2 = getCPUTime() - start;

	try {
		int bound = (int)ceil(model.get(GRB_DoubleAttr_ObjBound));

	}
	catch (GRBException e) {
		cout << "time limit in LB calculation" << endl;
		//std::cerr << "Gurobi Exception: " << e.getMessage() << std::endl;
		//std::cerr << "Error code: " << e.getErrorCode() << std::endl;
		sol.vLB2 = 0;
		return 0;
	}
	catch (...) {
		//std::cerr << "An unexpected error occurred." << std::endl;
		sol.vLB2 = 0;
		return 0;
	}
	sol.vLB2 = (int)ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);
	return (int)ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);
}


int LB_PCC_CONF(const Instance& inst, Solution& sol) {
	double start = getCPUTime();
	// compute bounds
	//int sumArea = 0;
	int LB = 0;
	int UB = 0;
	for (int j = 0; j < inst.n; j++) {
		//sumArea += inst.items[j][0] * inst.items[j][1] * inst.items[j][2];
		UB += inst.items[j][1] * inst.items[j][2];
	}
	//LB = sumArea / inst.W;
	LB = max(sol.vLB1, sol.vLB2); // start from the best known LB.
	// create node matrix
	vector<vector<int> > matToId(inst.items.size() + 1, vector<int>(UB + 1, -1));
	vector<vector<int> > idToMat;
	int nodeCount = 0;

	// create item arcs
	vector<bool> isA(UB + 1, false); isA[0] = true;
	vector<vector<int> > arcs;
	int tail, head;
	for (int j = 0; j < inst.items.size(); j++) {
		for (int i = UB; i >= 0; i--) {
			if (isA[i]) {
				if (matToId[j][i] == -1) {
					matToId[j][i] = nodeCount;
					idToMat.push_back({ j,i });
					nodeCount++;
				}
				tail = matToId[j][i];
				for (int k = 0; k <= inst.items[j][2]; k++) {
					if (i + inst.items[j][1] * k <= UB) {
						if (matToId[j + 1][i + inst.items[j][1] * k] == -1) {
							matToId[j + 1][i + inst.items[j][1] * k] = nodeCount;
							idToMat.push_back({ j + 1,i + inst.items[j][1] * k });
							nodeCount++;
						}
						head = matToId[j + 1][i + inst.items[j][1] * k];
						arcs.push_back({ tail,head,j,k });
						isA[i + inst.items[j][1] * k] = true;
					}
					else
						break;
				}
			}
		}
	}

	// Create a Gurobi environment

	GRBEnv env = GRBEnv();

	for (;;) {
		//cout << "Try H = " << LB << endl;

		// create a model           
		GRBModel model = GRBModel(env);
		model.getEnv().set(GRB_IntParam_OutputFlag, 0); //silence the output

		GRBLinExpr obj = 0;

		// declaration of the variables for the model
		vector<GRBVar> x(arcs.size());

		// initizalization of the variables for the model
		for (int k = 0; k < arcs.size(); k++) {
			if (idToMat[arcs[k][1]][1] <= LB)
				x[k] = model.addVar(0, inst.W, 0, GRB_CONTINUOUS);
		}
		model.update();

		// create linear expressions
		vector<GRBLinExpr> assigned(inst.items.size(), 0);
		vector<GRBLinExpr> cIn(idToMat.size(), 0);
		vector<GRBLinExpr> cOut(idToMat.size(), 0);
		for (int k = 0; k < arcs.size(); k++) {
			if (idToMat[arcs[k][1]][1] <= LB) {
				assigned[arcs[k][2]] += x[k] * arcs[k][3];
				cIn[arcs[k][1]] += x[k];
				cOut[arcs[k][0]] += x[k];
			}
		}

		model.update();

		// create assignment constraints
		for (int j = 0; j < inst.items.size(); j++)
			model.addConstr(assigned[j] == inst.items[j][2] * inst.items[j][0]);

		// create flow conservation constraints 
		for (int i = 1; i < idToMat.size(); i++) {
			if (idToMat[i][1] <= LB)
				model.addConstr(cIn[i] >= cOut[i]);
		}

		// set the objective: minimize z
		model.setObjective(cOut[0], GRB_MINIMIZE);

		// change some settings
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_IntParam_Method, 2);
		model.getEnv().set(GRB_IntParam_MIPFocus, 1);
		model.getEnv().set(GRB_DoubleParam_TimeLimit, max(1 - (getCPUTime() - start), 0.01));
		model.getEnv().set(GRB_DoubleParam_Cutoff, inst.W + 0.0001);

		// find the optimal solution		
		model.optimize();

		// store the results in a Solution object

		// if a solution has been found
		if (model.get(GRB_IntAttr_SolCount) >= 1) {
			UB = LB;
			sol.tLB += getCPUTime() - start;
			sol.tLB3 = getCPUTime() - start;
			sol.vLB3 = LB;
			return LB;
		}
		else {
			if (model.get(GRB_IntAttr_Status) == 9)
				break;
			else
				LB++;
		}
	}
	sol.tLB += getCPUTime() - start;
	sol.tLB3 = getCPUTime() - start;
	sol.vLB3 = LB;
	return LB;
}





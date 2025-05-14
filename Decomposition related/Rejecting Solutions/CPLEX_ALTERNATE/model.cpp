#include "model.h"
#include "IncumbentCallbackSPP.h"
#include <ilcplex/ilocplex.h>
using namespace std;
//BM REJECT DB ONE BY ONE (INSTEAD OF FAKE REJECT) 

void ALT(Instance& inst, Solution& sol) {
	double start = getCPUTime(); sol.startTime = start;
	int seed = 42; mt19937 rng(seed);

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
	if (n2 == 0) { sol.opt = 1; return; }
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

	int H = sol.iLB; // start from the last found LB (if it is the first time, then it is already eqaul to sol.iLB)

	bool keepSearching = true;
	while (keepSearching && getCPUTime() - start < inst.timeLimit) {

		IloEnv env;
		try {
			IloModel model(env);
			vector<vector<IloIntVar>> x;
			x.resize(inst.W, vector<IloIntVar>(n2));

			// initizalization of the variables for the model
			for (int j = 0; j < n2; j++) {
				for (int i = 0; i <= inst.W - items2[j][0]; i++) {
					if (NPs[j][i]) {
						isWActive[i] = true;
						x[i][j] = IloBoolVar(env);
					}
				}
			}

			// Linear expressions for assignment and height constraints
			IloExprArray assigned(env, n2);
			IloExprArray height(env, inst.W);

			for (int i = 0; i < inst.W; i++) height[i] = IloExpr(env);

			for (int j = 0; j < n2; j++) {
				assigned[j] = IloExpr(env);
				for (int i = 0; i <= inst.W - items2[j][0]; i++) {
					if (NPs[j][i]) {
						assigned[j] += x[i][j];
						for (int ip = i; ip < i + items2[j][0]; ip++) {
							height[ip] += items2[j][1] * x[i][j];
						}
					}
				}
			}

			// Assignment constraints
			for (int j = 0; j < n2; j++) {
				model.add(assigned[j] == 1).setName(("Assignment_" + to_string(j)).c_str());
				assigned[j].end();
			}

			// Height constraints
			for (int i = 0; i < inst.W; i++) {
				if (isWActive[i]) model.add(height[i] <= H).setName(("Height" + to_string(i)).c_str());
			}

			IloCplex cplex(model);
			cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0);
			cplex.setParam(IloCplex::Param::Threads, 1);
			cplex.setParam(IloCplex::Param::TimeLimit, max(inst.timeLimit - (getCPUTime() - start), EPSILON));
			cplex.setParam(IloCplex::Param::RandomSeed, 1);

			int desired_size = 100000;
			cplex.setParam(IloCplex::Param::MIP::Pool::Capacity, desired_size); // how many to store
			cplex.setParam(IloCplex::Param::MIP::Limits::Populate, desired_size); // how many to generate before stopping

			// Register the custom callback
			IncumbentCallbackSPP* lazyCallback = new IncumbentCallbackSPP(env, inst, sol, x, H, NPs);
			cplex.use(lazyCallback);

			// Optimize
			cplex.populate(); //cplex.solve();

			cout << "Solution status: " << cplex.getStatus() << endl;
			if (sol.opt == 1 || cplex.getStatus() == IloAlgorithm::Optimal) {
				cout << "Optimal solution found." << endl;
				sol.opt = 1;
				keepSearching = false; break;
			}
			else if (cplex.getStatus() == IloAlgorithm::Infeasible) {
				if (sol.NfailSl > 0) {
					cout << "There are MP solutions you skipped..." << endl;
					sol.opt = -1; break;
				}
				cout << "The problem is infeasible." << endl;
				H++; sol.LB += 1; sol.tLBinc = getCPUTime() - start;
				cout << "No feasible solution at this height. Height is increased to " << H + sol.prepacked_height << endl;
			}
			else if ((cplex.getStatus() == IloAlgorithm::Unknown) && (getCPUTime() - start +1 > inst.timeLimit) && (sol.opt ==0)){
				cout << "The problem reached time limit." << endl;
				break;
			}
			else {
				cout << "Other status: " << cplex.getStatus() << endl;
				break;
			}
		}
		catch (const IloException& e) {
			std::cerr << "Error: " << e << std::endl;
			sol.opt = -1;
			env.end();
			break;
		}
		catch (...) {
			std::cerr << "Unknown exception caught." << std::endl;
			env.end();
			break;
		}
		env.end();

	} //END WHILE

}

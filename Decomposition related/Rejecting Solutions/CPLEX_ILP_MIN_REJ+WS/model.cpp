#include "model.h"
#include "IncumbentCallbackSPP.h"
#include <ilcplex/ilocplex.h>
using namespace std;

void ILP_CPLEX_MIN_REJ(Instance& inst, Solution& sol) {
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

		// Define objective to minimize height H
		IloIntVar z(env, 0, IloIntMax, "z");

		// Height constraints
		for (int i = 0; i < inst.W; i++) {
			if (isWActive[i]) model.add(height[i] <= z).setName(("Height" + to_string(i)).c_str());
		}

		model.add(IloMinimize(env, z));

		IloCplex cplex(model);

		// ---- WARM START ----
		IloNumVarArray startVars(env);
		IloNumArray startVals(env);
		for (int i = 0; i < n2; i++) {
			int p_warm = get<0>(sol.packing_locations[i]);
			int q = get<1>(sol.packing_locations[i]);
			startVars.add(x[p_warm][i]);
			startVals.add(1); // 0 or 1 from your heuristic
		}
		cplex.addMIPStart(startVars, startVals);

		cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0);
		cplex.setParam(IloCplex::Param::Threads, 1);
		cplex.setParam(IloCplex::Param::TimeLimit, max(inst.timeLimit - (getCPUTime() - start), EPSILON));
		cplex.setParam(IloCplex::Param::RandomSeed, 1);
		cplex.setParam(IloCplex::Param::MIP::Limits::LowerObjStop, sol.iLB + 0.01); //Sets a lower objective stop. In a minimization MILP or MIQP, the solver will abort the optimization process as soon it finds a solution of value lower than or equal to the specified value for the new parameter lower objective stop. 

		//cplex.setParam(IloCplex::Param::MIP::Tolerances::UpperCutoff, sol.UB);
		//IloCplex::Param::MIP::Tolerances::UpperCutoff
		//Sets the upper cutoff tolerance.When the problem is a minimization problem, 
		//CPLEX cuts off or discards any solutions that are greater than the specified upper cutoff value.
		//If the model has no solution with an objective value less than or equal to the cutoff value, 
		//CPLEX declares the model infeasible.In other words, setting an upper cutoff value c for a minimization problem 
		//is similar to adding this constraint to the objective function of the model : obj <= c.

		// Register the custom callback
		IncumbentCallbackSPP lazyCallback(env, inst, sol, x, z, NPs);
		cplex.use(&lazyCallback);

		// Optimize
		cplex.solve();

		double bestLowerBound = ceil(cplex.getBestObjValue() - EPSILON); // When a model has been solved to optimality, this value matches the optimal solution value. Before optimality has been proven, this value is computed for a minimization (maximization) problem as the minimum (maximum) objective function value of all remaining unexplored nodes.
		//For a regular MIP optimization, this value is also the best known bound on the optimal solution value of the MIP problem. In fact, when a problem has been solved to optimality, this value matches the optimal solution value.
		sol.LB = bestLowerBound;
		cout << "Solution status: " << cplex.getStatus() << endl;
		if (cplex.getStatus() == IloAlgorithm::Optimal || sol.opt) {
			cout << "Optimal solution found." << endl;
			sol.opt = 1;

			// Retrieve the optimal objective value (minimized height)
			if (cplex.getStatus() == IloAlgorithm::Optimal) {
				double optimalHeight = cplex.getObjValue();
				cout << "Optimal objective value (minimized height): " << optimalHeight << endl;
				sol.UB = optimalHeight + sol.prepacked_height;
			}
			else sol.UB += sol.prepacked_height;
		}
		else if (cplex.getStatus() == IloAlgorithm::Infeasible) {
			cout << "The problem is infeasible." << endl;
		}
		else if (cplex.getStatus() == IloAlgorithm::Unknown) {
			cout << "The problem reached time limit." << endl;
		}
		else {
			cout << "Other status: " << cplex.getStatus() << endl;
			sol.UB += sol.prepacked_height;
		}

	}
	catch (const IloException& e) {
		std::cerr << "Error: " << e << std::endl;
		sol.opt = -1;
		env.end();
	}
	catch (...) {
		std::cerr << "Unknown exception caught." << std::endl;
		env.end();
	}
	env.end();

}

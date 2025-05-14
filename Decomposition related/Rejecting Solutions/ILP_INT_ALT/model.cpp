#include "model.h"


void ILP_INT_ALT(Instance& inst, Solution& sol) {
	double start = getCPUTime(); sol.startTime = start;
	int seed = 42; mt19937 rng(seed);
	// local variable
	vector<bool> isWActive(inst.W, false);

	// compute normal patterms
	vector<vector<bool> > NPs(inst.n, vector<bool>(inst.W, false));
	for (int j = 0; j < inst.n; j++) {
		if (inst.items[j][2] == 0) continue;
		NPs[j][0] = true;
		for (int jp = 0; jp < inst.n; jp++) {
			int lim = inst.items[jp][2];
			if (jp == j) lim--;
			for (int l = 0; l < lim; l++) {
				for (int i = inst.W - inst.items[jp][0] - 1; i >= 0; i--) {
					if (NPs[j][i]) NPs[j][i + inst.items[jp][0]] = true;
				}
			}
		}
	}
	inst.NPs = NPs;

	int H = sol.iLB;
	bool keepSearching = true;
	while (keepSearching && getCPUTime() - start + EPSILON < inst.timeLimit) {
		// create a model
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 1);
		try {
			GRBModel model = GRBModel(env);
			// declaration of the variables for the model
			vector<vector<GRBVar>> x;
			x.resize(inst.W, vector<GRBVar>(inst.n));

			// initizalization of the variables for the model
			for (int j = 0; j < inst.n; j++) {
				if (inst.items[j][2] == 0) continue;
				for (int i = 0; i <= inst.W - inst.items[j][0]; i++) {
					if (NPs[j][i]) {
						isWActive[i] = true;
						x[i][j] = model.addVar(0, inst.items[j][2], 0, GRB_INTEGER);
					}
				}
			}

			// create linear expressions
			vector<GRBLinExpr> assigned(inst.n, 0);
			vector<GRBLinExpr> height(inst.W, 0);
			for (int j = 0; j < inst.n; j++) {
				if (inst.items[j][2] == 0) continue;
				for (int i = 0; i <= inst.W - inst.items[j][0]; i++) {
					if (NPs[j][i]) {
						assigned[j] += x[i][j];
						for (int ip = i; ip < i + inst.items[j][0]; ip++) {
							height[ip] += inst.items[j][1] * x[i][j];
						}
					}
				}
			}
			// create assignment constraints
			for (int j = 0; j < inst.n; j++) {
				if (inst.items[j][2] == 0) continue;
				model.addConstr(assigned[j] == inst.items[j][2]);
			}

			// create height constraints 
			for (int i = 0; i < inst.W; i++) {
				if (isWActive[i])
					model.addConstr(height[i] <= H);
			}

			// change some settings
			model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
			model.getEnv().set(GRB_IntParam_Threads, 1);
			model.getEnv().set(GRB_DoubleParam_TimeLimit, max(inst.timeLimit - (getCPUTime() - start), EPSILON)); //use the remaining time left if you increase the LB
			model.getEnv().set(GRB_IntParam_Seed, 1);

			//////SETTINGS FOR GETTING MULTIPLE ALTERNATIVE SOLUTIONS
			// Limit how many solutions to collect
			model.set(GRB_IntParam_PoolSolutions, GRB_MAXINT); //defined to match the largest value that can be stored in an integer on your system, like 2147483647 (which is 2^31 - 1, the largest value of a signed 32-bit integer).
			// Limit the search space by setting a gap for the worst possible solution that will be accepted
			model.set(GRB_DoubleParam_PoolGap, 0);
			// do a systematic search for the k-best solutions
			model.set(GRB_IntParam_PoolSearchMode, 2);

			// create a callback
			mycallbackT2 cb = mycallbackT2(inst, sol, x, H, NPs);
			model.set(GRB_IntParam_LazyConstraints, 1); 		// indicate that we want to add Lazy Constraints
			model.setCallback(&cb);                     		// link the callback to the model

			// optimize	
			model.optimize();
			cout << "Model status = " << model.get(GRB_IntAttr_Status) << endl; //9,11,2

			if (sol.opt == 1) { // if a solution has been found (AFTER ADDING CUTS/ALTERNATING) note: it does not count a solution if it is not proved to be feasible in the callback.
				//aborted because opt is proven.
				cout << "Optimal found" << endl;
				keepSearching = false; break;
			}
			else if (sol.opt == 0 && sol.timeT < inst.timeLimit - 2 && model.get(GRB_IntAttr_Status) == 3) {
				if (sol.NfailSl > 0) { //this means you ignored solutions just because slave could not prove infeasibility, do not increase H. exit
					cout << "There are MP solutions you skipped. You need to come back to those before increasing H!" << endl;
					sol.opt = -1; keepSearching = false; break;
				}
				//otherwise, it means no solution of this height can be found.
				sol.tLBinc = getCPUTime() - start;
				H++; sol.LB += 1;
				cout << "No feasible solution at this height.\n Height is increased to = " << H << endl;
			}
			else {
				cout << "Time limit reached" << endl;
				keepSearching = false; break;
			}
		}
		catch (GRBException e) {
			std::cout << "Error code = " << e.getErrorCode() << std::endl;
			std::cout << e.getMessage() << std::endl;
			sol.opt = -1;
			break;
		}
		catch (...) {
			std::cout << "Exception during optimization" << std::endl;
			break;
		}

	}
}


void mycallbackT2::callback() {
	if (where == GRB_CB_MIPSOL) { // if a new incumbent is found
		if (getDoubleInfo(GRB_CB_RUNTIME) + 2.0 > inst.timeLimit) {
			cout << "Time limit reached: abort..." << endl;
			abort(); return;
		}
		if (sol.tInitSol == -1) sol.tInitSol = getCPUTime() - sol.startTime;

		sol.Ncalls += 1; // increase the number of calls by 1
		cout << "Ncalls " << sol.Ncalls << " calls " << endl;

		// transform the solution into a CSP solution by renumbering identical items (save the actual item id (number))
		vector<vector<int> > items_slave; //CSP
		int n_slave = 0;
		for (int j = 0; j < inst.n; j++) {
			if (inst.items[j][2] == 0) continue; //if it is already prepacked, skip
			for (int i = 0; i <= inst.W - inst.items[j][0]; i++) {
				if (NPs[j][i] && ceil(getSolution(x[i][j]) - EPSILON) > 0) {
					//cout << "(" << i << "," << j << ")= "<< ceil(getSolution(x[i][j]) - EPSILON) << endl;
					for (int l = 0; l < ceil(getSolution(x[i][j]) - EPSILON); l++) {
						//cout << "\t (" << i << "," << j << ") " << "Task " << j << " x-coord " << i << endl;
						items_slave.push_back({ inst.items[j][0],inst.items[j][1],i,j }); //item widht, height, position, org_id
						n_slave += 1;
					}
				}
			}
		}
		// store the solution of the master
		vector<vector<int>> bins(inst.W); //create a vector which stores the vertical item pieces packed in each column (bin)
		// get bin for each item
		for (int jp = 0; jp < n_slave; jp++) {
			//push back to an amount which is equal to its width
			int org_id = items_slave[jp][3]; //original item id
			int start_pos = items_slave[jp][2]; //starting position of the item acc to master solution
			for (int v = start_pos; v < start_pos + items_slave[jp][0]; v++) {
				bins[v].push_back(jp); //notice, you push back the item index in the model not the actual item index j (in other words not org_id)
			}
		}
		//print_master_sol(bins); // print the master solution per bin

		// CALL THE SLAVE
		vector<vector<int>> current_sol;
		double startSlave = getCPUTime(); //to measure the time spent in the slave
		int status = slave(items_slave, bins, H, current_sol, 0);
		if (status == 1) { // feasible 
			cout << "Slave feasible." << endl;
			sol.opt = 1; abort();
		}
		else {
			// SKIP THE SOLUTION
			cout << "Skipping the solution due to:";
			if (status == 0) { //time limit reached
				cout << "Time limit in slave..." << endl;
				sol.NfailSl += 1; //increase the nb fail in slave
			}
			if (status == 3) { //infeasible
				cout << "Slave infeasible..." << endl;
				sol.Ncuts++; // in this case this means number of solutions proven to be infeasible
			}
			cout << "Ncalls = " << sol.Ncalls << ": infeasible vs failed = " << sol.Ncuts << ", " << sol.NfailSl << endl;
			cout << "-------------------------------------------------------------------" << endl;
		}
		// END SLAVE //
		sol.tSlave += getCPUTime() - startSlave; //to measure the time spent in the slave
		cout << "Time Slave=" << getCPUTime() - startSlave << endl;
		cout << "-------------------------------------------------------------------" << endl;
		// Terminate Gurobi when time limit in callback is reached.
		cout << "Time info = " << getDoubleInfo(GRB_CB_RUNTIME) << endl;
		if (getDoubleInfo(GRB_CB_RUNTIME) + 2.0 > inst.timeLimit) {
			cout << "time limit reached: abort..." << endl;
			abort(); return;
		}
	}
}







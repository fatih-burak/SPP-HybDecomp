#include "model.h"


void ILP_BIN_ALT(Instance& inst, Solution& sol) {
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

	int H = sol.iLB;

	bool keepSearching = true;
	while (keepSearching && getCPUTime() - start < inst.timeLimit - 2) {
		// create a model
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 1);
		try {
			GRBModel model = GRBModel(env);
			// declaration of the variables for the model
			vector<vector<GRBVar> > x;
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
			model.update();

			// create linear expressions
			vector<GRBLinExpr> assigned(n2, 0);
			vector<GRBLinExpr> height(inst.W, 0);
			for (int j = 0; j < n2; j++) {
				for (int i = 0; i <= inst.W - items2[j][0]; i++) {
					if (NPs[j][i]) {
						assigned[j] += x[i][j];
						for (int ip = i; ip < i + items2[j][0]; ip++) {
							height[ip] += items2[j][1] * x[i][j];
						}
					}
				}
			}
			model.update();

			// create assignment constraints
			for (int j = 0; j < n2; j++) {
				model.addConstr(assigned[j] == 1);
			}

			vector<GRBConstr> c_height(inst.W);
			// create height constraints 
			for (int i = 0; i < inst.W; i++) {
				if (isWActive[i]) {
					c_height[i] = model.addConstr(height[i] <= H);
				}
			}
			model.update();

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
			mycallbackT1 cb = mycallbackT1(inst, sol, x, H, NPs);
			model.set(GRB_IntParam_LazyConstraints, 1); 		// indicate that we want to add Lazy Constraints
			model.setCallback(&cb);                     		// link the callback to the model

			// optimize	
			model.optimize();
			cout << model.get(GRB_IntAttr_Status) << endl; //9,11,2

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
		}
	}
	//calculate average similarity
	if (!sol.percentages.empty()) {
		double sum = 0.0;
		for (double num : sol.percentages) sum += num;
		sol.perc_sim = sum / sol.percentages.size();  // Calculate average
	}
	cout << "Percentage similarity = " << sol.perc_sim << endl;

}


void mycallbackT1::callback() {
	if (where == GRB_CB_MIPSOL) { // if a new incumbent is found
		if (getDoubleInfo(GRB_CB_RUNTIME) + 2.0 > inst.timeLimit) {
			cout << "Time limit reached: abort..." << endl;
			abort(); 
		}
		if (sol.tInitSol == -1) sol.tInitSol = getCPUTime() - sol.startTime;

		double start = getCPUTime(); //to measure the time spent in the slave
		sol.Ncalls += 1; // increase the number of calls by 1
		cout << "Ncalls " << sol.Ncalls << " calls " << endl;

		int n2 = inst.items2.size(); //nb of items
		vector<int> getStart(n2);
		for (int j = 0; j < n2; j++) {
			for (int i = 0; i <= inst.W - inst.items2[j][0]; i++) {
				if (NPs[j][i] && ceil(getSolution(x[i][j]) - EPSILON) > 0) getStart[j] = i;
			}
		}

		vector<vector<int> > items_slave(n2);
		vector<int> this_sol(n2); //! for measure similarity btw solutions
		for (int j = 0; j < n2; j++) {
			items_slave[j] = { inst.items2[j][0],inst.items2[j][1], getStart[j] }; //width, height, x-coord
			this_sol[j] = getStart[j]; //!
		}


		//!
		if (!sol.all_sols.empty()) {
			//this means there are solutions found before. lets compare the previous solution to this new solution
			//access the last one
			const vector<int>& last_vector = sol.all_sols.back();
			int count_similar = 0;
			for (int j = 0; j < n2; j++) {
				if (last_vector[j] == this_sol[j]) count_similar++;
			}
			sol.percentages.push_back(count_similar / (double)n2);
			cout << "Percent similar = " << count_similar / (double)n2 << endl;
		}

		//add this solution
		sol.all_sols.push_back(this_sol); //!

		vector<vector<int>> bins(inst.W); //create a vector which stores the vertical item pieces packed in each column (bin)
		for (int j = 0; j < n2; j++) {
			for (int k = getStart[j]; k < getStart[j] + items_slave[j][0]; k++) bins[k].push_back(j);
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
			abort();
		}
	}
}

#include "model.h"

void CP_ALT(Instance& inst, Solution& sol) {
	int start = getCPUTime();
	int seed = 42; mt19937 rng(seed);
	uniform_int_distribution<> dis(1, 10000);  // Range: [1, 10000]
	// from CSP to BPP
	vector<vector<int>> items2;
	int n2 = 0;
	for (int j = 0; j < inst.n; j++) {
		for (int k = 0; k < inst.items[j][2]; k++) {
			items2.push_back({ inst.items[j][0],inst.items[j][1],j });
			cout << "\titem=" << n2 << ": " << inst.items[j][0] << ", " << inst.items[j][1] << ", " << j << endl;
			n2++;
		}
	}
	inst.items2 = items2;
	int H = sol.iLB;

	// create a model
	IloEnv env;
	IloModel model(env);
	// declaration of the variables for the model
	IloIntervalVarArray x(env, n2);
	IloCumulFunctionExpr z(env);
	// initizalization of the variables for the model
	for (int j = 0; j < n2; j++) {
		x[j] = IloIntervalVar(env);
		x[j].setStartMin(0);
		x[j].setStartMax(inst.W - items2[j][0]);
		x[j].setSizeMin(items2[j][0]);
		x[j].setSizeMax(items2[j][0]);
		z += IloPulse(x[j], items2[j][1]);
	}

	// set the height constraint
	vector<IloConstraint> height_ILO_consts;
	IloConstraint height_constraint = (z <= H);
	height_ILO_consts.push_back(height_constraint);
	model.add(height_ILO_consts[0]);

	// change some settings
	IloCP cp(model);
	cp.setParameter(IloCP::Workers, 1);
	cp.setParameter(IloCP::LogPeriod, 500000);
	cp.setParameter(IloCP::FailureDirectedSearch, IloCP::On);
	cp.setParameter(IloCP::CumulFunctionInferenceLevel, IloCP::Extended);

	bool keepSearching = true; //CHECK
	while (keepSearching && getCPUTime() - start < inst.timeLimit - 1) {
		//SOLVE
		cp.setParameter(IloCP::SearchType, IloCP::Restart);
		cp.setParameter(IloCP::TimeLimit, max(inst.timeLimit - (getCPUTime() - start), 0.001));
		cp.setParameter(IloCP::RandomSeed, 1);
		//cp.setParameter(IloCP::LogVerbosity, IloCP::Quiet);
		cp.propagate();

		cp.startNewSearch(); // then in a while loop cp.next() to check all solutions
		while (cp.next() && getCPUTime() - start < inst.timeLimit-1) { //-2 seconds to ensure no funny business
			sol.Ncalls++;
			int master_status = cp.getStatus();
			cout << "\tNcalls = " << sol.Ncalls << endl;
			if (master_status == 1) { // a master solution is found
				if (sol.tInitSol == -1.0) sol.tInitSol = getCPUTime() - start;  // time spent until the very first master solution

				//get the current solution with relevant item info
				vector<vector<int> > items_slave(n2);
				vector<int> this_sol(n2); //! for measure similarity btw solutions
				for (int j = 0; j < n2; j++) {
					items_slave[j] = { inst.items2[j][0],inst.items2[j][1],(int)cp.getStart(x[j]) }; //width, height, x-coord
					this_sol[j] = (int)cp.getStart(x[j]); //!
					//cout << "Task-" << j << "=" << (int)cp.getStart(x[j]) << endl;
				}; 

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
					for (int k = cp.getStart(x[j]); k < cp.getEnd(x[j]); k++) bins[k].push_back(j);
					this_sol[j] = (int)cp.getStart(x[j]); //!
				}
				//print_master_sol(bins);

				// CALL THE SLAVE
				double startS = getCPUTime(); //to measure the time spent in in the slave
				vector<vector<int>> current_sol;
				int status = slave(items_slave, bins, H, current_sol, 0);
				cout << fixed << setprecision(2) << "Slave status = " << status << " and time in slave = " << getCPUTime() - startS << endl;
				if (status == 1) {
					cout << "\tSlave feasible." << endl;
					sol.solution = current_sol;
					cout << "\tOptimal found" << endl; sol.opt = 1;
					if (sol.tLBinc < 0) sol.tLBinc = 0;
					keepSearching = false; cp.endSearch(); //cp.end(); env.end();
					break;
				}
				else if (status == 0) { //time limit reached
					cout << "\tTime limit in slave: skipping the solution..." << endl;
					sol.NfailSl++; //increase the nb fail in slave (This number is not included in Ncuts)
					cout << "\tNfailSl added " << sol.NfailSl << " cuts..." << endl;
					cout << "\t-------------------------------------------------------------------" << endl;
				}
				else if (status == 3) { //infeasible
					sol.Ncuts++; //in this context this means number of proven infeasible MPs.
					cout << "\tSlave infeasible." << endl;
				}
				sol.tSlave += getCPUTime() - startS; //to measure the time spent in the slave (all callback)
			}
		}

		if (keepSearching && !(cp.next()) && getCPUTime() - start < inst.timeLimit-1) { //-1 seconds to ensure no funny business
			cout << cp.next() << " vs " << !(cp.next()) <<endl;
			if (sol.NfailSl > 0) {
				cout << "There are slave problems you skipped and before proving infeasibility you cannot increase the H!" << endl;
				sol.opt = -1; 
				keepSearching = false; cp.endSearch(); break;
			}
		 	if(getCPUTime() - start < inst.timeLimit -2){
				H++; //if infeasibility is proven and no MP solution is skipped and time limit is not exceeded
				sol.LB += 1; sol.tLBinc = getCPUTime() - start; //update the time to increase the LB
				cout << "No feasible solution at this height.\n Height is increased to = " << H << endl;
				model.remove(height_ILO_consts[0]); height_ILO_consts.clear();
				IloConstraint height_constraint = (z <= H); height_ILO_consts.push_back(height_constraint);
				model.add(height_ILO_consts[0]);
			}
			else break;
		}
		else break;
	}
	//calculate average similarity
	if (!sol.percentages.empty()) {
		double sum = 0.0;
		for (double num : sol.percentages) sum += num;
		sol.perc_sim = sum / sol.percentages.size();  // Calculate average
	}
	cout << "Percentage similarity = " << sol.perc_sim << endl;
}



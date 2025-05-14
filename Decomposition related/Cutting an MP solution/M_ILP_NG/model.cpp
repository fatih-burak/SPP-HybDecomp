#include "model.h"


void ILP_NG(Instance& inst, Solution& sol) {
	double start = getCPUTime();
	int seed = 42; mt19937 rng(seed); uniform_int_distribution<> dis(1, 10000);  // Range: [1, 10000]

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

	int H = sol.LB; // start from the last found LB (if it is the first time, then it is already eqaul to sol.iLB)
	cout << "Trying H = " << H + sol.prepacked_height << " actually trying H = " << H << endl; //show the tried height together with the prepacked height

	// create a model
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 1);
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

	// find the optimal solution  
	vector<GRBConstr> saved_GRBcuts;
	vector<vector<vector<int>>> forbid_start; //to store the current solution

	bool keepSearching = true;
	while (keepSearching && getCPUTime() - start < inst.timeLimit) {
		// add no-good cuts
		for (int cut = 0; cut < forbid_start.size(); cut++) {
			GRBLinExpr cut_this;
			int RHS = 0;
			int j = 0;
			while (1) {
				int item_no = forbid_start[cut][j][0];
				int previous_item = item_no;
				while (previous_item == item_no) {
					int item_start = forbid_start[cut][j][1];
					if (inst.NPs[item_no][item_start])	cut_this += x[item_start][item_no]; //add the index in the cut
					previous_item = item_no;
					j++;
					if (j == forbid_start[cut].size()) break;
					item_no = forbid_start[cut][j][0]; //move to the next item
				}
				RHS++;
				if (j == forbid_start[cut].size()) break;
			}
			saved_GRBcuts.push_back({ model.addConstr(cut_this <= RHS - 1) });
		}
		model.update();
		forbid_start.clear(); //empty the cuts list

		// change some settings
		model.getEnv().set(GRB_DoubleParam_TimeLimit, max(inst.timeLimit - (getCPUTime() - start), EPSILON)); //use the remaining time left
		//model.getEnv().set(GRB_IntParam_Seed, dis(rng)); // Random seed
		model.getEnv().set(GRB_IntParam_Seed, 1); // Fixed seed
		cout << "Current seed value: " << model.getEnv().get(GRB_IntParam_Seed) << std::endl;
		
		// optimize	
		model.optimize();

		// if a solution has been found
		if (model.get(GRB_IntAttr_Status) == 2) {// a master solution is found
			sol.Ncalls += 1; // increase the number of calls by 1
			cout << "Ncalls = " << sol.Ncalls << endl;
			if (sol.tInitSol == -1.0) sol.tInitSol = getCPUTime() - start;  // time spent until the very first master solution

			//get the current solution with relevant item info
			vector<int> getStart(n2);
			vector<int> this_sol(n2);//!
			for (int j = 0; j < n2; j++) {
				for (int i = 0; i <= inst.W - items2[j][0]; i++) {
					if (NPs[j][i] && x[i][j].get(GRB_DoubleAttr_X)) {
						getStart[j] = i;
						this_sol[j] = (int)i; //!
						//cout << "x[" << i << "][" << j << "] = 1;" << endl;
					}
				}
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

			vector<vector<int> > items_slave(n2);
			for (int j = 0; j < n2; j++) items_slave[j] = { inst.items2[j][0],inst.items2[j][1], getStart[j] }; //width, height, x-coord

			vector<vector<int>> bins(inst.W); //create a vector which stores the vertical item pieces packed in each column (bin)
			for (int j = 0; j < n2; j++) {
				for (int k = getStart[j]; k < getStart[j] + items_slave[j][0]; k++) bins[k].push_back(j);
			}
			//print_master_sol(bins);

			// CALL THE SLAVE
			double startS = getCPUTime(); //to measure the time spent in in the slave
			vector<vector<int>> current_sol;
			int status = slave(items_slave, bins, H, current_sol, 0);
			cout << fixed << setprecision(2) << "Slave status = " << status << " and time in slave = " << getCPUTime() - startS << endl;

			if (status == 1) {
				cout << "Slave feasible, optimal found!" << endl;
				sol.opt = 1;
				if (sol.tLBinc < 0) sol.tLBinc = 0;
				keepSearching = false;
			}
			else if (status == 0) { //time limit reached
				cout << "Time limit in slave: skipping the solution..." << endl;
				vector<vector<int>> cut; //to store the current solution
				for (int j = 0; j < n2; j++) cut.push_back({ j, getStart[j] }); // add the full cut
				forbid_start.push_back(cut);
				sol.NfailSl++; //increase the nb fail in slave (This number is not included in Ncuts)
				cout << "NfailSl added " << sol.NfailSl << " cuts..." << endl;
				cout << "-------------------------------------------------------------------" << endl;
			}
			else if (status == 3) { //infeasible
				cout << "Slave infeasible, adding the full cut (no-good cut)." << endl;

				vector<vector<int>> cut; //to store the current solution
				for (int j = 0; j < n2; j++) cut.push_back({ j, getStart[j] }); // add the full cut
				forbid_start.push_back(cut);
				sol.Ncuts++;

				cout << "Ncuts added " << sol.Ncuts << " cuts " << endl;
				cout << "-------------------------------------------------------------------" << endl;
			}
			sol.tSlave += getCPUTime() - startS; //to measure the time spent in the slave
		}
		else if (model.get(GRB_IntAttr_Status) == 9) {
			cout << "Time limit is reached." << endl; keepSearching = false;
		}
		else if (model.get(GRB_IntAttr_Status) == 3) { //Infeasible
			if (sol.NfailSl > 0) {
				cout << "This means you ignored solutions just because slave could not prove infeasibility, do not increase H. EXITING..." << endl;
				sol.opt = -1; break;
				//Note: so far never happened...
			}
			//this means no solution of this height can be found.
			H++; sol.LB = H; // keep H updated
			cout << "Trying H = " << H + sol.prepacked_height << " actually trying H = " << H << endl; //show the tried height together with the prepacked height

			forbid_start.clear(); //empty the cuts list
			//DELETE previous cuts from the model
			for (int cut = 0; cut < saved_GRBcuts.size(); cut++) model.remove(saved_GRBcuts[cut]);
			saved_GRBcuts.clear();
			//update the height constraints
			for (int i = 0; i < inst.W; i++) {
				if (isWActive[i]) c_height[i].set(GRB_DoubleAttr_RHS, H);
			}
			model.update();
			// Check the prefound cuts and pre-add them if they are still valid with this new height
			if (sol.cuts.size() > 0) {
				re_add_cuts(items2, inst.W, H, sol);
				for (set<vector<int>> prev_cut : sol.cuts) {
					vector<bool> deleted(n2, true);
					vector<vector<int> > items_slave(n2);
					//Add the cut
					vector<vector<int>> cut; //to store the current solution
					for (vector<int> i : prev_cut) cut.push_back({ i[1], i[0] });
					forbid_start.push_back(cut);
				}
			}
			sol.tLBinc = getCPUTime() - start;
		}
		else cout << "Program should not enter here." << endl;
	}
	//calculate average similarity
	if (!sol.percentages.empty()) {
		double sum = 0.0;
		for (double num : sol.percentages) sum += num;
		sol.perc_sim = sum / sol.percentages.size();  // Calculate average
	}
	cout << "Percentage similarity = " << sol.perc_sim << endl;
}


#include "model.h"


void CPM(Instance& inst, Solution& sol) {
	int start = getCPUTime();
	int seed = 42; mt19937 rng(seed);
	uniform_int_distribution<> dis(1, 10000);  // Range: [1, 10000]
	// from CSP to BPP
	vector<vector<int>> items2;
	int n2 = 0;
	for (int j = 0; j < inst.n; j++) {
		for (int k = 0; k < inst.items[j][2]; k++) {
			items2.push_back({ inst.items[j][0],inst.items[j][1],j });
			//cout << "item=" << n2 << ": " << inst.items[j][0] << ", " << inst.items[j][1] << ", " << j << endl;
			n2++;
		}
	}
	inst.items2 = items2;

	int H = sol.iLB; // start from the initial LB
	sol.LB = H; // starting H
	cout << "Trying H = " << H + sol.prepacked_height << " actually trying H = " << H << endl; //show the tried height together with the prepacked height

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
	cp.setParameter(IloCP::LogPeriod, 100000);
	cp.setParameter(IloCP::FailureDirectedSearch, IloCP::On);
	cp.setParameter(IloCP::CumulFunctionInferenceLevel, IloCP::Extended);

	// find the optimal solution  
	vector<IloConstraint> saved_ILOcuts;
	vector<vector<vector<int>>> forbid_start; //to store the current solution

	bool keepSearching = true; //CHECK
	while (keepSearching && getCPUTime() - start < inst.timeLimit) {
		if (sol.Ncalls >= 10) break;

		// Adding cuts
		for (int cut = 0; cut < forbid_start.size(); cut++) {

			//for (vector<int> cut_itself : forbid_start[cut]) {
			//	cout << "(" << cut_itself[0] << "," << cut_itself[1] << ")";
			//}
			//cout << " " << endl;

			IloOr OR_const(env);
			int j = 0;
			while (1) {
				IloAnd AND_const(env);
				int item_no = forbid_start[cut][j][0];
				int previous_item = item_no;
				while (previous_item == item_no) {
					//cout << "(" << item_no << "," << forbid_start[cut][j][1] << ")";
					int item_start = forbid_start[cut][j][1];
					AND_const.add(IloStartOf(x[item_no]) != item_start);
					previous_item = item_no;
					j++;
					if (j == forbid_start[cut].size()) break;
					item_no = forbid_start[cut][j][0]; //move to the next item
				}
				//cout << " " << endl;
				OR_const.add(AND_const);
				if (j == forbid_start[cut].size()) break;
			}
			model.add(OR_const);
			saved_ILOcuts.push_back(OR_const);
		}
		forbid_start.clear(); //empty the cuts list

		//SOLVE
		cp.setParameter(IloCP::SearchType, IloCP::Restart);
		cp.setParameter(IloCP::TimeLimit, max(inst.timeLimit - (getCPUTime() - start), 0.01));
		//cp.setParameter(IloCP::DefaultInferenceLevel, IloCP::Extended);
		// Generate a single random integer to set the seed
		int random_number = dis(rng);
		cp.setParameter(IloCP::RandomSeed, random_number);
		cp.setParameter(IloCP::LogVerbosity, IloCP::Quiet); 
		cp.propagate();
		cp.solve();

		int master_status = cp.getStatus();

		if (master_status == 1) { // a master solution is found
			sol.Ncalls += 1; // increase the number of calls by 1
			cout << "Ncalls = " << sol.Ncalls << endl;
			if (sol.tInitSol == -1.0) sol.tInitSol = getCPUTime() - start;  // time spent until the very first master solution

			//get the current solution with relevant item info
			vector<vector<int> > items_slave(n2);
			for (int j = 0; j < n2; j++) items_slave[j] = { inst.items2[j][0],inst.items2[j][1],(int)cp.getStart(x[j]) }; //width, height, x-coord

			vector<vector<int>> bins(inst.W); //create a vector which stores the vertical item pieces packed in each column (bin)
			for (int j = 0; j < n2; j++) {
				for (int k = cp.getStart(x[j]); k < cp.getEnd(x[j]); k++) bins[k].push_back(j);
			}
			//print_master_sol(bins);

			// CALL THE SLAVE
			double startS = getCPUTime(); //to measure the time spent in in the slave
			vector<vector<int>> current_sol;
			int status = slave(items_slave, bins, H, current_sol, 0);
			sol.tRealSl += getCPUTime() - startS; //to measure the time spent in the real slave
			cout << fixed << setprecision(2) << "Slave status = " << status << " and time in slave = " << getCPUTime() - startS << endl;

			//// BEGIN: SLAVE TEST
			double timeSpentSlave = getCPUTime();
			int status_CP = slave_CP(items_slave, bins, H, current_sol, 0, inst.name);
			int status_ILP = slave_ILP(items_slave, bins, H, current_sol, 0, inst.name);

			if ((status_CP == 1 && status_ILP == 3) || (status_ILP == 1 && status_CP == 3)) {
				for (int mim = 0; mim <= 10000; mim++) cout << inst.name << endl;
			}

			double addBackTime = getCPUTime() - timeSpentSlave;
			inst.timeLimit += addBackTime;
			cout <<"(Org CP, CP, ILP) = (" << status << ", " << status_CP << ", " << status_ILP << ")" << endl;

			//// END: SLAVE TEST

			if (status == 1) {
				cout << "Slave feasible." << endl;
				sol.solution = current_sol;
				sol.tFeasSl = getCPUTime() - startS; //update the stat

				cout << "Optimal found" << endl; sol.opt = 1;
				if (sol.tLBinc < 0) sol.tLBinc = 0;

				cout << "Container dimensions = " << inst.W << "x" << H << endl;
				//for (int j = 0; j < n2; j++) cout << "Rectangle-" << j << " of size=" << items_slave[j][0] << "x" << items_slave[j][1] << " x coordinate = " << current_sol[j][1] << ",y coordinate = " << current_sol[j][2] << endl;

				keepSearching = false;
				cp.end(); env.end();
			}
			else if (status == 0) { //time limit reached
				cout << "Time limit in slave: skipping the solution..." << endl;
				vector<vector<int>> cut; //to store the current solution
				for (int j = 0; j < n2; j++) cut.push_back({ j, (int)cp.getStart(x[j]) }); // add the full cut
				forbid_start.push_back(cut);
				sol.NfailSl++; //increase the nb fail in slave (This number is not included in Ncuts)
				cout << "NfailSl added " << sol.NfailSl << " cuts..." << endl;
				cout << "-------------------------------------------------------------------" << endl;
			}
			else if (status == 3) { //infeasible

				vector<bool> deleted(n2, false);
				int size_deleted = 0;

				cout << "Slave infeasible." << endl;
				double startMIS = getCPUTime(); //to measure the time spent in MIS finding
				bool atLeast1Cut = false; //to whether you add at least one cut. Otherwise, add the full cut  by lifting it

				// MIS-LR+O //
				double SmisLR = getCPUTime();
				deleted = MIS_LR(items_slave, bins, H);
				size_deleted = count(deleted.begin(), deleted.end(), true);
				if (size_deleted) { //only perform this operation if you have deleted something with MIS-LR, otherwise you will do this in the next step anyway!

					int sizeb4 = size_deleted;
					//perform MIS-O on top of MIS-LR
					deleted = MIS_O(items_slave, bins, H, deleted);
					size_deleted = count(deleted.begin(), deleted.end(), true);
					bool isValidCut = cut_check(items_slave, deleted, sol);
					if (isValidCut) {
						if (size_deleted - sizeb4 > 0) { sol.NmisLR_o++; sol.DmisLR_o += size_deleted; }
						vector<vector<int>> indices_MIS = LIFTING(inst, items_slave, deleted);
						//Add the cut
						vector<vector<int>> cut; //to store the current solution
						for (vector<int> i : indices_MIS) cut.push_back({ i[1], i[0] });
						forbid_start.push_back(cut);  atLeast1Cut = true;
						sol.Ncuts++; sol.NmisLR++; sol.DmisLR += sizeb4;
					}
				}
				sol.tmisLR += getCPUTime() - SmisLR;

				// MIS-O //
				double SmisO = getCPUTime();
				fill(deleted.begin(), deleted.end(), false); //reset it to all false (nothing is deleted)
				deleted = MIS_O(items_slave, bins, H, deleted);
				size_deleted = count(deleted.begin(), deleted.end(), true);
				if (size_deleted) { //if you do not delete anything, then do not enter since it is basically adding a full cut.
					sol.DmisO += size_deleted;
					bool isValidCut = cut_check(items_slave, deleted, sol);
					if (isValidCut) {
						vector<vector<int>> indices_MIS = LIFTING(inst, items_slave, deleted);
						//Add the cut
						vector<vector<int>> cut; //to store the current solution
						for (vector<int> i : indices_MIS) cut.push_back({ i[1], i[0] });
						forbid_start.push_back(cut); atLeast1Cut = true;
						sol.Ncuts++; sol.NmisO++;
					}
				}
				sol.tmisO += getCPUTime() - SmisO;

				// MIS-R //
				double SmisR = getCPUTime();
				fill(deleted.begin(), deleted.end(), false); //reset it to all false (nothing is deleted)
				int seed_value = dis(rng);
				deleted = MIS_R(items_slave, bins, H, deleted, seed_value);
				size_deleted = count(deleted.begin(), deleted.end(), true);
				if (size_deleted) { //if you do not delete anything, then do not enter since it is basically adding a full cut.
					sol.DmisR += size_deleted;
					bool isValidCut = cut_check(items_slave, deleted, sol);
					if (isValidCut) {
						vector<vector<int>> indices_MIS = LIFTING(inst, items_slave, deleted);
						//Add the cut
						vector<vector<int>> cut; //to store the current solution
						for (vector<int> i : indices_MIS) cut.push_back({ i[1], i[0] }); //(j,p)
						forbid_start.push_back(cut); atLeast1Cut = true;
						sol.Ncuts++; sol.NmisR++;
					}
				}
				sol.tmisR += getCPUTime() - SmisR;

				//If no MIS worked, add the full cut by lifting it
				if (!atLeast1Cut) {
					vector<vector<int>> indices_MIS = LIFTING(inst, items_slave, deleted);
					//Add the cut
					vector<vector<int>> cut; //to store the current solution
					for (vector<int> i : indices_MIS) cut.push_back({ i[1], i[0] });
					forbid_start.push_back(cut);
					sol.Ncuts++; sol.NmisFull++;
				}


				sol.tmis += getCPUTime() - startMIS;
				cout << "Ncuts added " << sol.Ncuts << " cuts " << endl;
				cout << "-------------------------------------------------------------------" << endl;
			}
			sol.tSlave += getCPUTime() - startS; //to measure the time spent in the slave (all callback)
		}

		else if (master_status == 0) { //time limit in master reached
			cout << "Time limit is reached." << endl;
			keepSearching = false;
			cp.end(); env.end();
		}
		else if (master_status == 3) { //master is infeasible
			if (sol.NfailSl > 0) {
				cout << "This means you ignored solutions just because slave could not prove infeasibility, do not increase H. EXITING..." << endl;
				sol.opt = 1903;
				break;
				//Note: check whether this ever happens. If yes, then you should go back and check the ignored slave problems. But so far never happened...
			}
			//this means no solution of this height can be found.
			H++; sol.LB = H; // keep H updated
			cout << "Trying H = " << H + sol.prepacked_height << " actually trying H = " << H << endl; //show the tried height together with the prepacked height
			forbid_start.clear(); //empty the cuts list
			//DELETE previous cuts from the model aka finding a relaxed solution whereas the height is not correct. (mostly class10: 10_020_02 etc)
			for (int cut = 0; cut < saved_ILOcuts.size(); cut++) model.remove(saved_ILOcuts[cut]);

			model.remove(height_ILO_consts[0]);
			height_ILO_consts.clear();
			IloConstraint height_constraint = (z <= H);
			height_ILO_consts.push_back(height_constraint);
			model.add(height_ILO_consts[0]);

			// Check the prefound cuts and pre-add them if they are still valid with this new height
			if (sol.cuts.size() > 0) re_add_cuts(items2, inst.W, H, sol);

			for (set<vector<int>> prev_cut : sol.cuts) {
				vector<bool> deleted(n2, true);
				vector<vector<int> > items_slave(n2);
				for (vector<int> i : prev_cut) {
					items_slave[i[0]] = { inst.items2[i[0]][0],inst.items2[i[0]][1], i[1] }; //width, height, x-coord
					deleted[i[0]] = false;
				}
				vector<vector<int>> indices_MIS = LIFTING(inst, items_slave, deleted);
				//Add the cut
				vector<vector<int>> cut; //to store the current solution
				for (vector<int> i : indices_MIS) cut.push_back({ i[1], i[0] });
				forbid_start.push_back(cut);
			}
			sol.tLBinc = getCPUTime() - start;
		}
		else cout << "Program should not enter here." << endl;
	}

}



#include "model.h"
void MP_CP(Instance& inst, Solution& sol) {
	double start = getCPUTime();
	random_device rd;  mt19937 gen(rd()); uniform_int_distribution<int> distrib(1, 10000); // Range [1, 10000] // Seed the generator
	mt19937 rng(distrib(gen)); uniform_int_distribution<> dis(1, 10000);  // Range: [1, 10000]
	
	// from CSP to BPP
	vector<vector<int>> items2;
	int n2 = 0;
	for (int j = 0; j < inst.n; j++) {
		for (int k = 0; k < inst.items[j][2]; k++) {
			items2.push_back({ inst.items[j][0],inst.items[j][1],j });
			//cout << "\titem=" << n2 << ": " << inst.items[j][0] << ", " << inst.items[j][1] << ", " << j << endl;
			n2++;
		}
	}
	inst.items2 = items2;

	int H = sol.LB; // start from the last found LB (if it is the first time, then it is already eqaul to sol.iLB)
	//cout << "Trying H = " << H + sol.prepacked_height << " actually trying H = " << H << endl; //show the tried height together with the prepacked height

	// create a model
	IloEnv env;
	IloModel model(env);
	// declaration of the variables for the model
	IloIntervalVarArray x(env, n2);
	IloCumulFunctionExpr z(env);
	// initizalization of the variables for the model
	for (int j = 0; j < n2; j++) {
		x[j] = IloIntervalVar(env);
		x[j].setName(("Task_" + to_string(j)).c_str());
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

	// find the optimal solution  
	vector<IloConstraint> saved_ILOcuts;
	vector<vector<vector<int>>> forbid_start; //to store the current solution

	//ADD THE CUTS COMING FROM PREVIOUS STATES
	if (sol.saved_cuts.size() > 0) cout << "Adding the cuts coming from previous states..." << endl;
	int numb_cut_added = 0;
	for (vector<vector<int>> prev_cut : sol.saved_cuts) {
		forbid_start.push_back(prev_cut);
		numb_cut_added++;
		cout << "Nb added = " << numb_cut_added << endl;
	}

	bool keepSearching = true; //CHECK
	while (keepSearching && getCPUTime() - start < inst.timeLimit) {
		// Adding cuts
		for (int cut = 0; cut < forbid_start.size(); cut++) {
			IloOr OR_const(env);
			int j = 0;
			while (1) {
				IloAnd AND_const(env);
				int item_no = forbid_start[cut][j][0];
				int previous_item = item_no;
				while (previous_item == item_no) {
					//cout << "\t(" << item_no << "," << forbid_start[cut][j][1] << ")";
					int item_start = forbid_start[cut][j][1];
					AND_const.add(IloStartOf(x[item_no]) != item_start);
					previous_item = item_no;
					j++;
					if (j == forbid_start[cut].size()) break;
					item_no = forbid_start[cut][j][0]; //move to the next item
				}
				//cout << "\t " << endl;
				OR_const.add(AND_const);
				if (j == forbid_start[cut].size()) break;
			}
			model.add(OR_const);
			saved_ILOcuts.push_back(OR_const);
		}
		forbid_start.clear(); //empty the cuts list

		//SOLVE
		cp.setParameter(IloCP::SearchType, IloCP::Restart);
		cp.setParameter(IloCP::TimeLimit, max(inst.timeLimit - (getCPUTime() - start), 0.0001));
		cp.setParameter(IloCP::LogVerbosity, IloCP::Quiet);
		int seedVal = dis(rng);
		cp.setParameter(IloCP::RandomSeed, seedVal); // Generate a single random integer to set the seed
		cout << "Current seed value: " << seedVal << std::endl;
		cout << "Time left in this model: " << max(inst.timeLimit - (getCPUTime() - start), 0.0001) << std::endl;

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
			if (status == 1) {
				cout << "Slave feasible, optimal is found." << endl;
				sol.solution = current_sol;
				sol.opt = 1;
				if (sol.tLBinc < 0) sol.tLBinc = 0;
				keepSearching = false;
				cp.end(); env.end();
			}
			else if (status == 0) { //time limit reached
				cout << "\tTime limit in slave: skipping the solution..." << endl;
				vector<vector<int>> cut; //to store the current solution
				for (int j = 0; j < n2; j++) cut.push_back({ j, (int)cp.getStart(x[j]) }); // add the full cut  //(j,p)
				forbid_start.push_back(cut);  sol.saved_cuts.push_back(cut);
				set<vector<int>> convertedCut(cut.begin(), cut.end()); sol.cuts.insert(convertedCut); //without lifting version

				sol.NfailSl++; //increase the nb fail in slave (This number is not included in Ncuts)
				cout << "\tNfailSl added " << sol.NfailSl << " cuts..." << endl;
				cout << "\t-------------------------------------------------------------------" << endl;
			}
			else if (status == 3) { //infeasible
				vector<bool> deleted(n2, false);
				int size_deleted = 0;

				cout << "\tSlave infeasible." << endl;
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
						vector<vector<int>> MIS_cut;
						for (int item_j = 0; item_j < items_slave.size(); item_j++) {
							if (!deleted[item_j]) MIS_cut.push_back({ item_j , items_slave[item_j][2] });//if j is not a deleted item,push the nonlifted item index
						}
						set<vector<int>> convertedCut(MIS_cut.begin(), MIS_cut.end()); sol.cuts.insert(convertedCut); //without lifting version
						//Now call LIFTING
						vector<vector<int>> indices_MIS = LIFTING(inst, items_slave, deleted);
						//Add the cut
						vector<vector<int>> cut; //to store the current solution
						for (vector<int> i : indices_MIS) cut.push_back({ i[1], i[0] });  //(j,p)
						forbid_start.push_back(cut);  atLeast1Cut = true;  sol.saved_cuts.push_back(cut);
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
						vector<vector<int>> MIS_cut;
						for (int item_j = 0; item_j < items_slave.size(); item_j++) {
							if (!deleted[item_j]) MIS_cut.push_back({ item_j , items_slave[item_j][2] });//if j is not a deleted item,push the nonlifted item index
						}
						set<vector<int>> convertedCut(MIS_cut.begin(), MIS_cut.end()); sol.cuts.insert(convertedCut); //without lifting version
						//Now call LIFTING
						vector<vector<int>> indices_MIS = LIFTING(inst, items_slave, deleted);
						//Add the cut
						vector<vector<int>> cut; //to store the current solution
						for (vector<int> i : indices_MIS) cut.push_back({ i[1], i[0] });  //(j,p)
						forbid_start.push_back(cut); atLeast1Cut = true;  sol.saved_cuts.push_back(cut);
						sol.Ncuts++; sol.NmisO++;
					}
				}
				sol.tmisO += getCPUTime() - SmisO;

				// MIS-R //
				double SmisR = getCPUTime();
				fill(deleted.begin(), deleted.end(), false); //reset it to all false (nothing is deleted)
				//int seed_value = distrib(gen);
				int seed_value = dis(rng);
				deleted = MIS_R(items_slave, bins, H, deleted, seed_value);
				size_deleted = count(deleted.begin(), deleted.end(), true);
				if (size_deleted) { //if you do not delete anything, then do not enter since it is basically adding a full cut.
					sol.DmisR += size_deleted;
					bool isValidCut = cut_check(items_slave, deleted, sol);
					if (isValidCut) {
						vector<vector<int>> MIS_cut;
						for (int item_j = 0; item_j < items_slave.size(); item_j++) {
							if (!deleted[item_j]) MIS_cut.push_back({ item_j , items_slave[item_j][2] });//if j is not a deleted item,push the nonlifted item index
						}
						set<vector<int>> convertedCut(MIS_cut.begin(), MIS_cut.end()); sol.cuts.insert(convertedCut); //without lifting version
						//Now call LIFTING
						vector<vector<int>> indices_MIS = LIFTING(inst, items_slave, deleted);
						//Add the cut
						vector<vector<int>> cut; //to store the current solution
						for (vector<int> i : indices_MIS) cut.push_back({ i[1], i[0] }); //(j,p)
						forbid_start.push_back(cut); atLeast1Cut = true;  sol.saved_cuts.push_back(cut);
						sol.Ncuts++; sol.NmisR++;
					}
				}
				sol.tmisR += getCPUTime() - SmisR;

				//If no MIS worked, add the full cut by lifting it
				if (!atLeast1Cut) {
					vector<vector<int>> MIS_cut;
					for (int item_j = 0; item_j < items_slave.size(); item_j++) {
						if (!deleted[item_j]) MIS_cut.push_back({ item_j , items_slave[item_j][2] });//if j is not a deleted item,push the nonlifted item index
					}
					set<vector<int>> convertedCut(MIS_cut.begin(), MIS_cut.end()); sol.cuts.insert(convertedCut); //without lifting version
					//Now call LIFTING
					vector<vector<int>> indices_MIS = LIFTING(inst, items_slave, deleted);
					//Add the cut
					vector<vector<int>> cut; //to store the current solution
					for (vector<int> i : indices_MIS) cut.push_back({ i[1], i[0] });
					forbid_start.push_back(cut);  sol.saved_cuts.push_back(cut);
					sol.Ncuts++; sol.NmisFull++;
				}

				sol.tmis += getCPUTime() - startMIS;
				cout << "\tNcuts added " << sol.Ncuts << " cuts " << endl;
				cout << "\t-------------------------------------------------------------------" << endl;
			}
			sol.tSlave += getCPUTime() - startS; //to measure the time spent in the slave (all callback)
		}
		else if (master_status == 0) { //time limit in master reached
			cout << "\tTime limit is reached." << endl;
			keepSearching = false;
			cp.end(); env.end();
		}
		else if (master_status == 3) { //master is infeasible
			if (sol.NfailSl > 0) {
				cout << "This means you ignored solutions just because slave could not prove infeasibility, do not increase H. EXITING..." << endl;
				sol.opt = -1; break;
				//Note: so far never happened...
			}
			//this means no solution of this height can be found.
			cout << "\tincreasing H" << endl;
			H++; sol.LB = H; // keep H updated
			//cout << "Trying H = " << H + sol.prepacked_height << " actually trying H = " << H << endl; //show the tried height together with the prepacked height
			sol.tLBinc = getCPUTime() - start;
			sol.isLBincreased = true;

			forbid_start.clear(); //empty the cuts list
			sol.saved_cuts.clear(); //empty the grand cuts list

			//DELETE previous cuts from the model aka finding a relaxed solution whereas the height is not correct. (mostly class10: 10_020_02 etc)
			for (int cut = 0; cut < saved_ILOcuts.size(); cut++) model.remove(saved_ILOcuts[cut]);

			model.remove(height_ILO_consts[0]);
			height_ILO_consts.clear();
			IloConstraint height_constraint = (z <= H);
			height_ILO_consts.push_back(height_constraint);
			model.add(height_ILO_consts[0]);

			// Check the prefound cuts and pre-add them if they are still valid with this new height
			if (sol.cuts.size() >= 1) {
				re_add_cuts(items2, inst.W, H, sol);

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
					for (vector<int> i : indices_MIS) cut.push_back({ i[1], i[0] }); //(j,p)
					forbid_start.push_back(cut); sol.saved_cuts.push_back(cut);
				}
			}
			keepSearching = false; break; // EXIT and REWARD RL
		}
		else cout << "\tProgram should not enter here." << endl;
	}
}



#include "slave_functions.h"

void re_add_cuts(const vector<vector<int>>& items2, const int& W, const int& H, Solution& sol) {
	cout << "\tReadding the cuts by checking them first with the new height..." << endl;
	set<set<vector<int>>> new_cuts;
	set<set<vector<int>>> new_skipped;
	int count_success = 0;
	int count_skipped = 0;
	int size_of_cuts = sol.cuts.size();
	// Checks the prefound cuts in this new DB of height H
	for (set<vector<int>> cut : sol.cuts) {
		vector<vector<int> > items_slave(cut.size());
		int item_id = 0;
		for (vector<int> item : cut) {
			int j = item[0], p = item[1];
			//cout << "j=" << j << " p=" << p << endl;
			items_slave[item_id] = { items2[j][0], items2[j][1], p }; //width, height, x-coord
			item_id++;
		}
		vector<vector<int>> bins(W); //create a vector which stores the vertical item pieces packed in each column (bin)
		for (int jp = 0; jp < items_slave.size(); jp++) {
			for (int k = items_slave[jp][2]; k < items_slave[jp][2] + items_slave[jp][0]; k++) bins[k].push_back(jp);
		}

		double startS = getCPUTime(); //to measure the time spent in in the slave
		vector<vector<int>> current_sol;
		int status = slave(items_slave, bins, H, current_sol, 0);
		sol.tRealSl += getCPUTime() - startS; //to measure the time spent in the real slave
		//cout << fixed << setprecision(2) << "\t\t Prefound cut slave status = " << status << " and time in slave = " << getCPUTime() - startS << endl;

		if (status == 3) { //infeasible
			new_cuts.insert(cut); //(j,p)
			count_success++;
		}
		if (status == 0) { //time limit reached ones. save them directly to the saved_cuts because these will not be lifted!
			vector<vector<int>> modified_cut; 
			for (vector<int> item : cut) modified_cut.push_back({ item[0], item[1]}); //(j,p)
			sol.saved_cuts.push_back(modified_cut);
			count_skipped++;
		}

	}

	cout << "\tOut of " << size_of_cuts << " cuts, " << count_success << " cuts are still infeasible and " <<count_skipped << " are skipped!" << endl;
	sol.cuts = new_cuts;
	return;
}

bool cut_check(const vector<vector<int>>& items_slave, vector<bool>& deleted, Solution& sol) {
	double cutStart = getCPUTime();

	int n_slave = items_slave.size(); //number of items (n)
	set<vector<int>> new_cut;
	for (int j = 0; j < n_slave; j++) {
		if (deleted[j]) continue; //if j is a deleted item, continue
		new_cut.insert({ j,items_slave[j][2] });
	}

	/*
	* if checking super set is time-wise costly, just check whether you already found the "exact" cut and remove checking the bool statement
	auto it = sol.cuts.find(new_cut);
	if (it != sol.cuts.end()) {
		cout << "\tCut is already found!"<< endl;
		return false;
	}
	*/

	bool isSubset = false;
	for (set<vector<int>> cut : sol.cuts) {
		isSubset = includes(new_cut.begin(), new_cut.end(), cut.begin(), cut.end());
		if (isSubset) break;
	}
	if (isSubset) {
		cout << "\t\t- Cut is a super set of an already found cut, do not add!" << endl;
		sol.tCutCheck += getCPUTime() - cutStart;
		sol.NelCut++;
		return false;
	}
	else {
		sol.cuts.insert(new_cut);
	}
	sol.tCutCheck += getCPUTime() - cutStart;
	return true; //cut is not found before
}

vector<bool> MIS_LR(const vector<vector<int>>& items_slave, const vector<vector<int>>& bins, const int& H) {
	cout << "\tMIS-LR: left-right delete..." << endl;
	int n_slave = items_slave.size(); //number of items (n)
	int W = bins.size(); //size of the container

	//left delete
	vector<int> to_be_deleted; //empty list
	vector<bool> deleted_b(n_slave);

	vector<int> p_list;	vector<int> p_rev_list;
	for (int p = 0; p < bins.size(); p++) p_list.push_back(p);
	for (int p = bins.size() - 1; p >= 0; p--) p_rev_list.push_back(p);

	int p = 0; // the column
	int p_index = 0;

	bool reverse = false;
	for (;;) {
		if (reverse) p = p_rev_list[p_index];
		else p = p_list[p_index];

		if (bins[p].size() == 0) { p_index++; continue; }// nothing to delete

		int size1 = to_be_deleted.size(); // get the initial size of the vector
		for (int i : bins[p]) to_be_deleted.push_back(i); //this has to be an ordered list (smallest to largest)
		sort(to_be_deleted.begin(), to_be_deleted.end());
		to_be_deleted.erase(unique(to_be_deleted.begin(), to_be_deleted.end()), to_be_deleted.end());
		int size2 = to_be_deleted.size();

		if (size1 == size2) { p_index++;  continue; } //nothing new added

		vector<vector<int>> temp_items_slave;
		vector<vector<int>> temp_bins(bins.size()); // keep the reduced vertical positions of items
		int tracker = 0;
		int new_item_index = 0;
		// get bin for each item
		for (int jp = 0; jp < n_slave; jp++) {
			if (jp == to_be_deleted[tracker]) { //skip the ones that are going to be deleted
				tracker++;
				continue;
			}
			//push back to an amount which is equal to its width
			int start_pos = items_slave[jp][2]; //starting position of the item acc to master solution
			int width = items_slave[jp][0];
			for (int v = start_pos; v < start_pos + width; v++) {
				temp_bins[v].push_back(new_item_index); //notice, you push back the item index in the model not the actual item index j
			}
			temp_items_slave.push_back(items_slave[jp]);
			new_item_index++;
		}
		if (new_item_index == 0) { p_index++;  continue; }//nothing to check

		vector<vector<int>> dummy; //find a way to remove this...
		int status = slave(temp_items_slave, temp_bins, H, dummy, 1);

		if (status == 1 || status == 0) { //stop if either feasible is found or time limit is reached
			if (reverse) break;
			reverse = true;
			p_index = 0;
		}
		else {
			for (int i : to_be_deleted) deleted_b[i] = true;
			p_index++;
		}
	}

	int MIS1size = count(deleted_b.begin(), deleted_b.end(), true);
	cout << "\t\t Nb of items deleted in MIS-LR: " << MIS1size << endl;

	return deleted_b;
}

vector<bool> MIS_O(const vector<vector<int>>& items_slave, const vector<vector<int>>& bins, const int& H, vector<bool>& deleted) {
	cout << "\tMIS-O: ordered list (smallest to largest) delete..." << endl;
	int size_b4 = count(deleted.begin(), deleted.end(), true); //MIS size before
	int n_slave = items_slave.size(); //number of items (n)
	int W = bins.size(); //size of the container
	vector<vector<int>> remainingItems;
	for (int jp = 0; jp < n_slave; jp++) {
		if (!deleted[jp]) { //item is not deleted
			remainingItems.push_back({ items_slave[jp][0],items_slave[jp][1],jp });
		}
	}

	//sort the remaining items
	sort(remainingItems.begin(), remainingItems.end(), sortItemsArea);
	//for (auto i : remainingItems) cout << i[0] << " " << i[1] << " " << i[2] << endl;

	for (vector<int> item_info : remainingItems) {
		int j = item_info[2]; // assign the item index
		vector<vector<int>> temp_items_slave;
		vector<vector<int>> temp_bins(bins.size()); // keep the reduced vertical positions of items
		int new_item_index = 0;
		// get bin for each item
		for (int jp = 0; jp < n_slave; jp++) {
			if (jp == j || deleted[jp]) continue; //skip the one that is going to be tested and that are already deleted.

			//push back to an amount which is equal to its width
			int start_pos = items_slave[jp][2]; //starting position of the item acc to master solution
			int width = items_slave[jp][0];
			for (int v = start_pos; v < start_pos + width; v++) {
				temp_bins[v].push_back(new_item_index); //notice, you push back the item index in the model not the actual item index j
			}
			temp_items_slave.push_back(items_slave[jp]);
			new_item_index++;
		}
		if (new_item_index == 0) { continue; }//nothing to check

		vector<vector<int>> dummy; //find a way to remove this...
		int status = slave(temp_items_slave, temp_bins, H, dummy, 1);
		if (status != 1 && status != 0) deleted[j] = true;
		else break; //I think it is reasonable to end the loop if exclusion of a small item turns out to be feasible... (the larger one will be most probably feasible as well anyway?)
	}

	int size_af = count(deleted.begin(), deleted.end(), true); //MIS size after
	//cout << "\tdeleted items in MIS-2: "; for (auto o : deleted) cout << o << " "; cout << endl;
	cout << "\t\t Nb of items deleted in MIS-O = " << (size_af - size_b4) << endl;

	return deleted;
}

vector<bool> MIS_R(const vector<vector<int>>& items_slave, const vector<vector<int>>& bins, const int& H, vector<bool>& deleted, int seed_value) {
	cout << "\tMIS-R: random delete..." << endl;
	int size_b4 = count(deleted.begin(), deleted.end(), true); //MIS size before
	int n_slave = items_slave.size(); //number of items (n)
	int W = bins.size(); //size of the container

	vector<int> remainingItems;
	for (int jp = 0; jp < n_slave; jp++) {
		if (!deleted[jp]) remainingItems.push_back(jp);
	}

	mt19937 rng(seed_value);
	shuffle(remainingItems.begin(), remainingItems.end(), rng);
	//for (int i : remainingItems) cout << i << endl;

	double start_MIS_R = getCPUTime();

	for (int j : remainingItems) {
		vector<vector<int>> temp_items_slave;
		vector<vector<int>> temp_bins(bins.size()); // keep the reduced vertical positions of items
		int new_item_index = 0;
		// get bin for each item
		for (int jp = 0; jp < n_slave; jp++) {
			if (jp == j || deleted[jp]) continue; //skip the one that is going to be tested and that are already deleted.

			//push back to an amount which is equal to its width
			int start_pos = items_slave[jp][2]; //starting position of the item acc to master solution
			int width = items_slave[jp][0];
			for (int v = start_pos; v < start_pos + width; v++) {
				temp_bins[v].push_back(new_item_index); //notice, you push back the item index in the model not the actual item index j
			}
			temp_items_slave.push_back(items_slave[jp]);
			new_item_index++;
		}
		if (new_item_index == 0) { continue; }//nothing to check

		vector<vector<int>> dummy; //find a way to remove this...
		int status = slave(temp_items_slave, temp_bins, H, dummy, 1);

		if (status != 1 && status != 0) deleted[j] = true;

		if (getCPUTime() - start_MIS_R > 0.5) break;  //limit the number of times you will try to delete the items by limiting the time

	}
	int size_af = count(deleted.begin(), deleted.end(), true); //MIS size after
	//cout << "\tdeleted items in MIS-2: "; for (auto o : deleted) cout << o << " "; cout << endl;
	cout << "\t\t Nb of items deleted in MIS-R = " << (size_af - size_b4) << endl;
	return deleted;

}

vector<vector<int>> LIFTING(Instance& inst, const vector<vector<int>>& items_slave, vector<bool>& deleted) {
	int n_slave = items_slave.size(); //number of items (n)
	int W = inst.W; //size of the container
	int deleted_size = count(deleted.begin(), deleted.end(), true);

	// LIFTING

	vector<set<int>> vert_overlap(n_slave); //subset of items that vertically overlap with j K_s(j): I include j itself also in the list
	for (int j = 0; j < n_slave; j++) {
		if (deleted[j]) continue; //if j is a deleted item, continue
		int start_j = items_slave[j][2]; //starting position of the item acc to master solution
		int end_j = start_j + items_slave[j][0] - 1;
		for (int jp = 0; jp < n_slave; jp++) {
			if (j == jp || deleted[jp]) continue;
			int start_jp = items_slave[jp][2]; //starting position of the item acc to master solution
			int end_jp = start_jp + items_slave[jp][0] - 1;
			if ((start_jp >= start_j && start_jp <= end_j) || (start_j >= start_jp && start_j <= end_jp))  vert_overlap[j].insert(jp);
		}
	}

	// Create a Gurobi environment
	GRBEnv env = GRBEnv();

	env.set(GRB_IntParam_OutputFlag, 0);
	GRBModel model = GRBModel(env);
	// declaration of the variables for the model
	vector<GRBVar> l(n_slave); vector<GRBVar> r(n_slave);
	// initizalization of the variables for the model

	for (int j = 0; j < n_slave; j++) {
		if (deleted[j]) continue; //if jp is a deleted item, continue
		l[j] = model.addVar(0, items_slave[j][2], 0, GRB_CONTINUOUS); //lb=0, ub=p_j
		r[j] = model.addVar(items_slave[j][2], W - items_slave[j][0], 0, GRB_CONTINUOUS); //lb=p_j, ub=W-w_j
	}
	model.update();

	//create objective
	GRBLinExpr objective;
	for (int j = 0; j < n_slave; j++) {
		if (deleted[j]) continue; //if jp is a deleted item, continue
		objective += r[j] - l[j];
	}
	model.setObjective(objective, GRB_MAXIMIZE);
	model.update();

	// create constraints
	for (int j = 0; j < n_slave; j++) {
		if (deleted[j]) continue; //if jp is a deleted item, continue
		for (int i : vert_overlap[j]) {
			model.addConstr(l[j] + items_slave[j][0] >= r[i] + 1);
		}
	}
	model.update();

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_DoubleParam_TimeLimit, 1); //use 1 sec
	// optimize	
	model.optimize();

	vector<vector<int>> x_indices;
	for (int j = 0; j < n_slave; j++) {
		if (deleted[j]) continue; //if jp is a deleted item, continue

		int left = l[j].get(GRB_DoubleAttr_X); int right = r[j].get(GRB_DoubleAttr_X);
		//cout << "\tfor item " << j << ": l=" << left << " ,r=" << right << endl;
		for (int p = left; p <= right; p++) {
			x_indices.push_back({ p, j }); //x[p,j]
		}
		//To deactivate LIFTING and pushing only non-lifted item indices
		//cout << "\tfor item " << j << ": i=" << items_slave[j][2] << endl;
		//x_indices.push_back({ items_slave[j][2], j });
	}
	//cout << "\t\t Cut size is lifted to = " << x_indices.size() << endl;
	return x_indices;
}

int slave(const vector<vector<int>>& items_slave_org, const vector<vector<int>>& bins_org, const int& H, vector<vector<int>>& current_sol, int type) {
	//type: Preprocessing beginning of the algorithm = 2, MIS call=1, Main Slave call=0 , preprocess slave call = -1;

	vector<vector<int>> items_slave; vector<vector<int>> bins;
	items_slave = items_slave_org; bins = bins_org;
	if (type != -1) preprocess_slave2(items_slave, bins);

	int n_slave = items_slave.size();
	IloEnv env;
	IloModel model(env);
	// declaration of the variables for the model
	IloIntervalVarArray x(env, n_slave);
	// initizalization of the variables for the model
	for (int j = 0; j < n_slave; j++) {
		x[j] = IloIntervalVar(env);
		x[j].setStartMin(0);
		x[j].setStartMax(H - items_slave[j][1]); //min-max starting positions of the interval var 
		// extra constraint.
		//x[j].setEndMin(items2[j][0]); x[j].setEndMax(inst.W); //min-max ending positions of the interval var
		x[j].setSizeMin(items_slave[j][1]); x[j].setSizeMax(items_slave[j][1]); // size of the interval var
	}
	// no-overlap constraints
	IloArray<IloIntervalVarArray> items_in_column(env, bins.size());
	for (int i = 0; i < bins.size(); i++) {
		if (i > 0 && bins[i] == bins[i - 1]) continue; // Do not check the columns that are the same as one before
		items_in_column[i] = IloIntervalVarArray(env, bins[i].size());
		for (int j = 0; j < bins[i].size(); j++) {
			items_in_column[i][j] = x[bins[i][j]];
		}
		model.add(IloNoOverlap(env, items_in_column[i]));
	}
	// change some settings
	IloCP cp(model);
	if (type == 1) cp.setParameter(IloCP::TimeLimit, 0.5); //MIS call
	else if (type == 2) cp.setParameter(IloCP::TimeLimit, 0.3); //actual preprocess call (from the very beginning)
	else if (type == -1) cp.setParameter(IloCP::TimeLimit, 0.5); // slave preprocess call (maybe reduce, check the time)
	else cp.setParameter(IloCP::TimeLimit, 0.5); // actual slave

	cp.setParameter(IloCP::LogPeriod, 10000);
	cp.setParameter(IloCP::SearchType, IloCP::Restart);
	cp.setParameter(IloCP::Workers, 1);
	cp.setParameter(IloCP::LogVerbosity, IloCP::Quiet); //make slave silent
	cp.setParameter(IloCP::RandomSeed, 0);
	cp.solve();
	int status = cp.getStatus();

	if (status == 1) { // feasible
		// store the results 
		for (int j = 0; j < n_slave; j++) current_sol.push_back({ j, items_slave[j][2], (int)cp.getStart(x[j]) }); // (item id, x-coord, y-coord)

	}
	if (cp.getInfo(IloCP::SolveTime) > 0.25) {
		if (type == 1) cout << fixed << setprecision(2) << "\t\t mis slave = " << status << " and time = " << cp.getInfo(IloCP::SolveTime) << endl;
		if (type == -1) cout << fixed << setprecision(2) << "\t\t\t\t slave preprocess-1 = " << status << " and time = " << cp.getInfo(IloCP::SolveTime) << endl;
	}
	return status;
}

void preprocess_slave2(vector<vector<int>>& items_slave, vector<vector<int>>& bins) {
	int n_slave = items_slave.size();
	int W = bins.size();

	//stats
	int NSINc_w = 0; int NSDec_W = 0;

	//find the left right items
	vector<vector <int>> left(n_slave);
	vector<vector <int>> right(n_slave);
	for (int j = 0; j < n_slave; j++) {
		int start_j = items_slave[j][2];
		int width_j = items_slave[j][0];
		for (int i = 0; i < n_slave; i++) {
			if (j == i) continue;
			int start_i = items_slave[i][2];
			int width_i = items_slave[i][0];
			if (start_i + width_i <= start_j) left[j].push_back(i);
			if (start_i >= start_j + width_j) right[j].push_back(i);
		}
	}

	// increase item widths
	vector<vector<int>> remainingItems; //save items with their indices
	for (int jp = 0; jp < n_slave; jp++)
		remainingItems.push_back({ items_slave[jp][0],items_slave[jp][1], items_slave[jp][2],jp });

	//consider items one at a time by nondecreasing width (smallest to largest)
	sort(remainingItems.begin(), remainingItems.end(), sortWidthSL);

	vector <int> l_j(n_slave, 0); vector <int> r_j(n_slave, W); //left and right points of the items.
	for (vector<int> item : remainingItems) {
		int j = item[3];
		//find left most position that item j can be extended
		for (int i : left[j]) {
			int start_i = items_slave[i][2];
			int width_i = items_slave[i][0];
			if (l_j[j] < start_i + width_i) l_j[j] = start_i + width_i;
		}
		items_slave[j][2] = l_j[j];
		//find right most position that item j can be extended
		for (int i : right[j]) {
			int start_i = items_slave[i][2];
			if (r_j[j] > start_i) r_j[j] = start_i;
		}
		//if (items_slave[j][0] != r_j[j] - l_j[j]) { cout << "\t\t\t\t In the slave preprocess, item " << j << " is enlarged by " << r_j[j] - l_j[j] - items_slave[j][0] << endl; NSINc_w++; }
		items_slave[j][0] = r_j[j] - l_j[j];
	}

	//change bins!
	for (int i = 0; i < bins.size(); i++) bins[i].clear();
	for (int j = 0; j < n_slave; j++) {
		for (int k = items_slave[j][2]; k < items_slave[j][2] + items_slave[j][0]; k++) bins[k].push_back(j);
	}
}

//additional functions
void print_master_sol(const vector<vector<int>>& bins) {
	// print the solution per bin
	for (int i = 0; i < bins.size(); i++) {
		cout << "\tBin " << i << "\t [";
		for (int j = 0; j < bins[i].size(); j++)
			cout << bins[i][j] << " ";
		cout << "\t]" << endl;
	}
}



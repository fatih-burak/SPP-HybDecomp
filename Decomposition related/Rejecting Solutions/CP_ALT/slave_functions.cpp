#include "slave_functions.h"

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



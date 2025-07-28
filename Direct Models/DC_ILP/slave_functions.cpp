#include "slave_functions.h"

int slave(const vector<vector<int>>& items_slave_org, const vector<vector<int>>& bins_org, const int& H, vector<vector<int>>& current_sol, int type) {
	//type: Preprocessing beginning of the algorithm = 2, MIS call=1, Main Slave call=0 , preprocess slave call = -1;

	vector<vector<int>> items_slave; vector<vector<int>> bins;
	items_slave = items_slave_org; bins = bins_org;
	if (type != -1) preprocess_slave(items_slave, bins);

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
	if (type == 1) cp.setParameter(IloCP::TimeLimit, 1); //MIS call
	else if (type == 2) cp.setParameter(IloCP::TimeLimit, 0.3); //actual preprocess call (from the very beginning)
	else if (type == -1) cp.setParameter(IloCP::TimeLimit, 10); // slave preprocess call (maybe reduce, check the time)
	else cp.setParameter(IloCP::TimeLimit, 1); // actual slave

	cp.setParameter(IloCP::LogPeriod, 10000);
	cp.setParameter(IloCP::SearchType, IloCP::Restart);
	cp.setParameter(IloCP::Workers, 1);
	cp.setParameter(IloCP::LogVerbosity, IloCP::Quiet); //make slave silent
	cp.solve();
	int status = cp.getStatus();

	if (status == 1) { // feasible
		// store the results 
		for (int j = 0; j < n_slave; j++) current_sol.push_back({ j, items_slave[j][2], (int)cp.getStart(x[j]) }); // (item id, x-coord, y-coord)

	}
	return status;
}

void preprocess_slave(vector<vector<int>>& items_slave, vector<vector<int>>& bins) {
	int n_slave = items_slave.size();
	int W = bins.size();

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

	vector<vector<int>> remainingItems; //save items with their indices because you will sort them
	for (int jp = 0; jp < n_slave; jp++) remainingItems.push_back({ items_slave[jp][0],items_slave[jp][1], items_slave[jp][2],jp });

	// PREPROCESS SLAVE -1 :MERGE ITEMS
	//consider items one at a time by nonincreasing order of their starting position acc to master sol (largest to smallest)
	sort(remainingItems.begin(), remainingItems.end(), sortItemsStartPosLS);
	//for (auto item : remainingItems) cout << "Item-"<<item[3] << ": "<< item[0] << " " << item[1] << " " << item[2] << endl;

	vector<bool> isItemMerged(n_slave, false);
	for (vector<int> item : remainingItems) {
		//SELECT THE CANDIDATE BIG ITEM AND DO ALL THESE OPERATIONS FOR TO MERGE SOME SMALLER ITEMS WITH THIS ONE
		int j = item[3], w_j = item[0], h_j = item[1], p_j = item[2];
		if (isItemMerged[j]) continue; //if the item is already merged with another item, skip.

		///////////////////////////////////////////////////////////////////////////////////////
		//FIRST STEP IS TO CHECK WHETHER YOU CAN MERGE WITH ALL THE THINGS TO THE "LEFT" OF j//
		///////////////////////////////////////////////////////////////////////////////////////

		bool all_h_smaller = false; //check whether all h_i<=h_j for all i in left[j]
		for (int i : left[j]) {
			if (isItemMerged[i]) continue; //probably not necessary
			int h_i = items_slave[i][1];
			if (h_i > h_j) { all_h_smaller = false; break; }
			all_h_smaller = true;
		}
		bool isFit = false;
		// Check if all the items in left[j] fit into the substrip
		if (all_h_smaller) { //all h_i <= h_j for all i in left[j], attempt packing the items in left[j] into a substrip of width = p_j, height = h_j
			//CALL SLAVE FOR THE REDUCED STRIP
			int W_reduced_strip = p_j, H_reduced_strip = h_j;
			// create a sub_item list (of size left[j]) with relevant item info (by reindexing them)
			vector<vector<int> > sub_items;
			for (int i : left[j]) {
				if (isItemMerged[i]) continue;
				sub_items.push_back({ items_slave[i][0],items_slave[i][1],items_slave[i][2] }); //width, height, x-coord
			}
			if (sub_items.size() == 0) continue;

			vector<vector<int>> sub_bins(W_reduced_strip); //create a vector which stores the vertical item pieces packed in each column (bin)
			for (int jp = 0; jp < sub_items.size(); jp++) {
				for (int k = sub_items[jp][2]; k < sub_items[jp][2] + sub_items[jp][0]; k++) sub_bins[k].push_back(jp);
			}
			//print_master_sol(sub_bins);

			double startS = getCPUTime(); //to measure the time spent in in the slave
			vector<vector<int>> dummy;
			int status = slave(sub_items, sub_bins, H_reduced_strip, dummy, -1);
			//cout << fixed << setprecision(2) << "Slave preprocess : merge items status = " << status << " and time it took = " << getCPUTime() - startS << endl;
			if (status == 1) {
				//cout << "\tLeft: Managed to merge item " << j << "(larger), with items: ";
				//mark the items as merged
				//isItemMerged[j] = true; // do not mark the big item as merged since you will declare this as the new merged item and you want to check this again on your way to RIGHT merge
				for (int i : left[j]) {
					if (isItemMerged[i]) continue;
					isItemMerged[i] = true; //mark the smaller items as merged
					//cout << i << " "; 
				}
				isFit = true;
				// I think it is okay to modify the item-j and make it the larger item
				item = { p_j + w_j,h_j,0 }; //probably not ness
				items_slave[j] = { p_j + w_j,h_j,0 };
				//cout << endl;
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////
		//SECOND STEP IS TO CHECK WHETHER YOU CAN MERGE WITH ALL THE THINGS TO THE "RIGHT" OF j//
		/////////////////////////////////////////////////////////////////////////////////////////

		all_h_smaller = false; //check whether all h_i<=h_j for all i in right[j]
		for (int i : right[j]) {
			if (isItemMerged[i]) continue;
			int h_i = items_slave[i][1];
			if (h_i > h_j) { all_h_smaller = false; break; }
			all_h_smaller = true;
		}

		isFit = false;
		// Check if all the items in right[j] fit into the substrip
		if (all_h_smaller) { //all h_i <= h_j for all i in right[j], attempt packing the items in right[j] into a substrip of width = p_j, height = h_j
			//CALL SLAVE FOR THE REDUCED STRIP
			int W_reduced_strip = W - (p_j + w_j), H_reduced_strip = h_j;
			// create a sub_item list (of size right[j]) with relevant item info (by reindexing them)
			vector<vector<int> > sub_items;
			for (int i : right[j]) {
				if (isItemMerged[i]) continue;
				sub_items.push_back({ items_slave[i][0],items_slave[i][1],items_slave[i][2] - (p_j + w_j) }); //width, height, x-coord //re-adjust the coordinate 
			}

			if (sub_items.size() == 0) continue;

			vector<vector<int>> sub_bins(W_reduced_strip); //create a vector which stores the vertical item pieces packed in each column (bin)
			for (int jp = 0; jp < sub_items.size(); jp++) {
				for (int k = sub_items[jp][2]; k < sub_items[jp][2] + sub_items[jp][0]; k++) sub_bins[k].push_back(jp);
			}
			//print_master_sol(sub_bins);

			double startS = getCPUTime(); //to measure the time spent in in the slave
			vector<vector<int>> dummy;
			int status = slave(sub_items, sub_bins, H_reduced_strip, dummy, -1);
			//cout << fixed << setprecision(2) << "Slave preprocess : merge items status = " << status << " and time it took = " << getCPUTime() - startS << endl;
			if (status == 1) {
				//cout << "\tRight: Managed to merge item " << j << "(larger), with items: ";
				//mark the items as merged
				//isItemMerged[j] = true; // do not mark the big item as merged since you will declare this as the new merged item and you want to check this again on your way to RIGHT merge
				for (int i : right[j]) {
					if (isItemMerged[i]) continue;
					isItemMerged[i] = true;//mark the smaller items as merged
					//cout << i << " ";
				}
				isFit = true;
				item = { W - p_j, h_j, p_j }; //probably not ness
				items_slave[j] = { W - p_j, h_j, p_j };
				//cout << endl;
			}
		}

	}
	vector<vector<int>> newlyCreatedItems; //save merged items 
	for (int jp = 0; jp < n_slave; jp++) {
		if (isItemMerged[jp]) continue;
		newlyCreatedItems.push_back({ items_slave[jp][0],items_slave[jp][1], items_slave[jp][2],jp });
		//cout << "Item-" << jp << ": " << items_slave[jp][0] << " " << items_slave[jp][1] << " " << items_slave[jp][2] << endl;
	}
	items_slave = newlyCreatedItems; //change the list
	n_slave = items_slave.size();
	//change bins!
	bins.clear(); bins.resize(W);
	for (int j = 0; j < newlyCreatedItems.size(); j++) {
		for (int k = items_slave[j][2]; k < items_slave[j][2] + items_slave[j][0]; k++) bins[k].push_back(j);
	}


	// PREPROCESS SLAVE -2 : INCREASE ITEM WIDTHS
	//consider items one at a time by nondecreasing width (smallest to largest)
	n_slave = items_slave.size();
	//find the left right items
	left.clear(); left.resize(n_slave);
	right.clear(); right.resize(n_slave);
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
	remainingItems.clear();
	//remainingItems.resize(n_slave);
	for (int jp = 0; jp < n_slave; jp++)
		remainingItems.push_back({ items_slave[jp][0],items_slave[jp][1], items_slave[jp][2],jp });


	sort(remainingItems.begin(), remainingItems.end(), sortWidthSL);
	vector <int> l_j(n_slave, 0); vector <int> r_j(n_slave, W); //left and right points of the items.
	for (vector<int> item : remainingItems) {
		int j = item[3];
		//find left most position that item j can be extended
		for (int i : left[j]) {
			int start_i = items_slave[i][2]; //if you activate preprocess -1, this is no longer correct
			int width_i = items_slave[i][0];
			if (l_j[j] < start_i + width_i) l_j[j] = start_i + width_i;
		}
		items_slave[j][2] = l_j[j];
		//find right most position that item j can be extended
		for (int i : right[j]) {
			int start_i = items_slave[i][2];
			if (r_j[j] > start_i) r_j[j] = start_i;
		}
		//if (items_slave[j][0] != r_j[j] - l_j[j]) { cout << "\t\t\t In the slave preprocess, item " << j << " is enlarged by " << r_j[j] - l_j[j] - items_slave[j][0] << endl; NSINc_w++; }
		items_slave[j][0] = r_j[j] - l_j[j];
	}

	//change bins!
	for (int i = 0; i < bins.size(); i++) bins[i].clear();
	for (int j = 0; j < n_slave; j++) {
		for (int k = items_slave[j][2]; k < items_slave[j][2] + items_slave[j][0]; k++) bins[k].push_back(j);
	}

}


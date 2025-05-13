#include "helper_functions.h"
#include "slave_functions.h"

int getValueForInstance(const string& filename, const string& instance_name) {
	ifstream file(filename); // Open the file
	if (!file.is_open()) {
		cerr << "Error: Could not open the file!" << endl;
		return -1; // Return an error code
	}

	string line;
	while (getline(file, line)) { // Read the file line by line
		istringstream lineStream(line);
		string current_instance_name;
		int value;

		lineStream >> current_instance_name >> value; // Parse instance name and value
		if (current_instance_name == instance_name) {
			return value; // Return the value if instance names match
		}
	}

	file.close();
	cerr << "Error: Instance name not found in the file!" << endl;
	return -1; // Return an error code if instance name is not found
}


void Instance::print() {
	cout << "n = " << n << " " << "W = " << W << endl;
	for (int j = 0; j < n; j++)
		cout << "Item " << j << ": " << items[j][0] << " " << items[j][1] << " " << items[j][2] << endl;
}

double getCPUTime() {
	return (double)clock() / CLOCKS_PER_SEC;
}

bool sortItems(const vector<int>& v1, const vector<int>& v2) {
	return M * v1[0] + v1[1] > M * v2[0] + v2[1];
}

bool sortItemsHeight(const vector<int>& v1, const vector<int>& v2) { //largest to smallest
	return M * v1[1] + v1[0] > M * v2[1] + v2[0];
}


bool sortItemsArea(const vector<int>& v1, const vector<int>& v2) { //largest to smallest
	return v1[1] * v1[0] < v2[1] * v2[0];
}


// Custom comparison function for sorting by the first element (width)
bool sortWidthSL(const std::vector<int>& v1, const std::vector<int>& v2) { //smallest to largest
	if (v1[0] == v2[0]) { //if tie: choose the smallest height
		return v1[1] < v2[1];
	}
	return v1[0] < v2[0];
}

// Custom comparison function for sorting by the third element (packed position acc to master sol)
bool sortItemsStartPosLS(const std::vector<int>& v1, const std::vector<int>& v2) { //largest to smallest
	if (v1[2] == v2[2]) { //if tie: choose the largest height first
		return v1[1] > v2[1];
	}
	return v1[2] > v2[2];
}

// Custom comparison function for sorting by the third element (packed position acc to master sol)
bool sortItemsStartPosSL(const std::vector<int>& v1, const std::vector<int>& v2) { //smallest to largest
	if (v1[2] == v2[2]) { //if tie: choose the largest width first
		return v1[0] > v2[0];
	}
	return v1[2] < v2[2];
}



Instance readInstance(const string& filename, Solution& sol) {
	// define variables
	Instance inst;

	// open the file
	ifstream file(filename);

	// read the file
	if (file.is_open()) { //if the file is open
		string line;

		// first line contains number of items
		getline(file, line, '\n');
		inst.n = stoi(line);
		sol.assignedBin.resize(inst.n);

		// reshape the array that will hold the items
		inst.items.resize(inst.n, vector<int>(3, 0));

		// second line contains the number of bins
		getline(file, line, '\n');
		inst.W = stoi(line); inst.Wo = inst.W;

		// the remaining lines contain the width, the height, and the demand of the items
		for (int j = 0; j < inst.n; j++) {
			getline(file, line, ' ');
			inst.items[j][0] = stoi(line);
			getline(file, line, ' ');
			inst.items[j][1] = stoi(line);
			getline(file, line, '\n');
			inst.items[j][2] = stoi(line);
		}

		// close the file
		file.close();
	}
	else {
		// if the file cannot be opened: print error and return default Inst
		cout << "Unable to open file";
	}
	// Preprocessing 1 -- prepack items	
	sort(inst.items.begin(), inst.items.end(), sortItems);

	//inst.print();
	Instance inst2 = inst; //will keep track of which items we are testing
	int sumAreaL = 0; int sumAreaR = 0;
	int heightL = 0; int heightR = 0;
	int idx2 = inst.n;
	for (int j = 0; j < inst.n; j++) inst2.items[j][2] = 0; //set all demand to 0

	for (int j = 0; j < inst.n; j++) {
		if (inst.items[j][0] <= inst.W / 2.0) break; //if an item is smaller than W/2, break the whole loop since the next ones will be even smaller (ordered list).

		//take large item j
		// The left
		heightL += inst.items[j][1] * inst.items[j][2]; //pack all of its copies on top of each other
		sumAreaL += inst.items[j][0] * inst.items[j][1] * inst.items[j][2];
		inst2.items[j][2] = inst.items[j][2]; //set the demand of item-j to its original demand (keeping track of what is packed)

		// The right
		while (idx2 >= j + 2 && inst.items[idx2 - 1][0] + inst.items[j][0] <= inst.W) { //starting from the last item (smallest), go backwards.
			idx2--;
			heightR = max(heightR, inst.items[idx2][1]); //assuming you are putting all the demanded amount next to each other (but they can be put also on top of each other,which will be considered in the feasibility test)
			sumAreaR += inst.items[idx2][0] * inst.items[idx2][1] * inst.items[idx2][2];
			inst2.items[idx2][2] = inst.items[idx2][2];
		}

		// Feasibility test
		//cout << "Test " << j << " " << idx2;
		if (sumAreaL + sumAreaR <= inst.W * heightL && heightL >= heightR) { //if there is enough area to pack and height of the big items are larger than the small items
			if (CP(inst2, heightL, sol.assignedBin) == 0) {
				sol.prepacked_height += heightL;
				sumAreaL = 0; sumAreaR = 0; heightL = 0; heightR = 0;
				for (int k = 0; k < inst.n; k++) inst2.items[k][2] = 0; //empty the helper inst
				for (int k = 0; k <= j; k++) inst.items[k][2] = 0; //delete the large item
				for (int k = idx2; k < inst.n; k++) inst.items[k][2] = 0; //delete the 
				//cout << " prepacked!" << endl;
			}
			else {
				//cout << " not prepacked" << endl;
			}
		}
		else {
			//cout << " no hope" << endl;
		}
	}

	sort(inst.items.begin(), inst.items.end(), sortItems);

	// Preprocessing 2 -- decrease W
	vector<bool> isR(inst.W + 1, false); isR[0] = true;
	for (int j = 0; j < inst.n; j++) {
		for (int l = 0; l < inst.items[j][2]; l++) {
			for (int i = inst.W - inst.items[j][0]; i >= 0; i--) {
				if (isR[i]) isR[i + inst.items[j][0]] = true;
			}
		}
	}
	//cout << "W from " << inst.W << " to ";
	while (!isR[inst.W] && inst.W >= 2) inst.W--;
	//cout << inst.W << endl;

	// Preprocessing 3 -- inrease w
	for (int j = 0; j < inst.n; j++) {
		if (inst.items[j][2] == 0) continue;
		isR.resize(0); isR.resize(inst.W - inst.items[j][0] + 1, false); isR[0] = true;
		for (int jp = 0; jp < inst.n; jp++) {
			int lim = inst.items[jp][2];
			if (jp == j) lim--;
			for (int l = 0; l < lim; l++) {
				for (int i = inst.W - inst.items[j][0] - inst.items[jp][0]; i >= 0; i--) {
					if (isR[i]) isR[i + inst.items[jp][0]] = true;
				}
			}
		}
		int tw = inst.W - inst.items[j][0];
		//cout << "Item " << j << " from " << inst.items[j][0] << " to ";
		while (!isR[tw]) tw--;
		if (2 * inst.items[j][0] > inst.W) inst.items[j][0] += (inst.W - inst.items[j][0] - tw);
		else inst.items[j][0] += (inst.W - inst.items[j][0] - tw) / inst.items[j][2];
		//cout << inst.items[j][0] << endl;
	}
	//inst.print();
	return inst;
}

int CP(const Instance& inst, const int& heightL, vector<vector<int> >& assignedBin) {

	// create a model
	IloEnv env;
	IloModel model(env);

	// from CSP to BPP
	vector<vector<int> > items2;
	int n2 = 0;
	for (int j = 0; j < inst.n; j++) {
		for (int k = 0; k < inst.items[j][2]; k++) {
			items2.push_back({ inst.items[j][0],inst.items[j][1],j });
			n2++;
		}
	}

	// declaration of the variables for the model
	IloIntervalVarArray x(env, n2);
	IloCumulFunctionExpr z(env);

	// initizalization of the variables for the model
	for (int j = 0; j < n2; j++) {
		x[j] = IloIntervalVar(env);
		x[j].setStartMin(0);
		if (items2[j][0] > inst.W / 2) x[j].setStartMax(0); //large items start at 0
		else x[j].setStartMax(inst.W - items2[j][0]);
		if (j > 0 && items2[j][2] == items2[j - 1][2])  model.add(IloStartBeforeStart(env, x[j - 1], x[j])); //identical item ordering
		x[j].setSizeMin(items2[j][0]); //size of the interval is the width of the item
		x[j].setSizeMax(items2[j][0]);
		z += IloPulse(x[j], items2[j][1]); //pulse height
	}

	// set the objective: minimize z
	model.add(z <= heightL);

	// change some settings
	IloCP cp(model);
	cp.setParameter(IloCP::TimeLimit, 0.1);
	cp.setParameter(IloCP::Workers, 1);
	cp.setParameter(IloCP::LogPeriod, 500000);
	cp.setParameter(IloCP::LogVerbosity, IloCP::Quiet);
	cp.solve();

	if (cp.getStatus() == 1) {

		vector<vector<int>> bins(inst.W); //create a vector which stores the vertical item pieces packed in each column (bin)
		// get bin for each item
		for (int j = 0; j < n2; j++) {
			for (int k = cp.getStart(x[j]); k < cp.getEnd(x[j]); k++) {
				bins[k].push_back(j); //notice, you push back the item index in the model not the actual item index items2[j][2]
			}
		}

		//get the current solution with relevant item info
		vector<vector<int> > items_slave(n2);
		for (int j = 0; j < n2; j++) items_slave[j] = { items2[j][0], items2[j][1],(int)cp.getStart(x[j]) }; //width, height, x-coord

		//call the y-check
		vector<vector<int>> dummy;
		int statusNew = slave(items_slave, bins, heightL, dummy, 2); //y-check
		if (statusNew == 1) { //if it is feasible, save in the assigned bin vector (you do preprocessing in slave, so rethink how you will get back to the actual coordinates).
			/*
			for (int j = 0; j < n2; j++) {
				// cout << "Task "<< j << " starts at time " << cp.getStart(x[j]) << endl;
				//assignedBin[items2[j][2]].push_back(cp.getStart(x[j]));
			}
			*/
			return 0;
		}
		else return -1;
	}
	if (cp.getStatus() == 0) {
		cout << " time limit in master preprocessing, skipping..." << endl;
	}
	return -1; //master is infeasible (or time limit reached)
}

void printInfo(const string& pathAndFileout, const Solution& sol, const string& filein) {
	string nameFile = pathAndFileout;
	std::ofstream file(nameFile.c_str(), std::ios::out | std::ios::app);
	file << fixed << setprecision(2) << filein << "\t" << sol.opt << "\t" << sol.timeT << "\t"
		<< sol.iLB << "\t" << sol.LB << "\t" << sol.UB << "\t" << sol.contR << endl;
	file.close();
}

void printSInfo(const Solution& sol, const Instance& inst) {
	vector<vector<int> > bins(inst.Wo);
	vector<int>	load(inst.Wo, 0);
	for (int j = 0; j < sol.assignedBin.size(); j++) {
		for (int i = 0; i < sol.assignedBin[j].size(); i++) {
			for (int ip = sol.assignedBin[j][i]; ip < sol.assignedBin[j][i] + inst.items[j][0]; ip++) {
				load[ip] += inst.items[j][1];
				bins[ip].push_back(j);
			}
		}
	}
	for (int i = 0; i < inst.Wo; i++) {
		cout << "Bin " << i << "\t" << load[i] << "\t [";
		for (int j = 0; j < bins[i].size(); j++)
			cout << bins[i][j] << " ";
		cout << "]" << endl;
	}
}

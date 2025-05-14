#include "helper_functions.h"
#include "slave_functions.h"
#include <cmath>

void readPolicyVectorFromFile(Agent& agent, const int& policyNo, const string& pathPolicyFile){
	ifstream file(pathPolicyFile);
	if (!file.is_open()) {
		cerr << "Error opening file!" << std::endl;
		return;
	}
	vector<int> policyData;
	string line;
	while (getline(file, line)) {
		stringstream stream(line);
		int currentPolicyNo;
		stream >> currentPolicyNo;
		// Check if the current line is the one we're looking for
		if (currentPolicyNo == policyNo) {
			// Read the rest of the values in the line into the vector
			int value;
			while (stream >> value)	policyData.push_back(value);
			file.close();
			break;
		}
	}

	cout << "Policy-" << policyNo << endl;
	for (auto i : policyData) cout << i << " ";
	cout << endl;

	// create the empty ones
	int nbLevels = 3, nbState = 7, nbAct = 3;
	// Initialize the 3D vector with dimensions 5x7x2 filled with 0.0
	vector< vector< vector<double>>> Q(nbLevels, vector< vector<double>>(nbState, vector<double>(nbAct, 0.0)));
	// Modify Q[0] to have only 1 sub-vector of dimension 3
	Q[0] = vector< vector<double>>(1, vector<double>(nbAct, 0.0));
	Q[1] = vector< vector<double>>(nbState, vector<double>(nbAct, 0.0));
	Q[2] = vector< vector<double>>(nbState, vector<double>(2, 0.0));

	// Initialize the 3D vector with dimensions 5x7x2 filled with 0
	vector< vector< vector<int>>> N(nbLevels, vector< vector<int>>(nbState, vector<int>(nbAct, 0)));
	N[0] = vector< vector<int>>(1, vector<int>(nbAct, 0));
	N[1] = vector< vector<int>>(nbState, vector<int>(nbAct, 0));
	N[2] = vector< vector<int>>(nbState, vector<int>(2, 0));
	agent.Q = Q; agent.N = N;

	agent.Q[0][0][policyData[0]] = 1;

	agent.Q[1][0][policyData[1]] = 1;
	agent.Q[1][1][policyData[2]] = 1;
	agent.Q[1][2][policyData[3]] = 1;
	agent.Q[1][3][policyData[4]] = 1;
	agent.Q[1][4][policyData[5]] = 1;
	agent.Q[1][5][policyData[6]] = 1;
	agent.Q[1][6][policyData[7]] = 1;

	agent.Q[2][0][policyData[8]] = 1;
	agent.Q[2][1][policyData[9]] = 1;
	agent.Q[2][2][policyData[10]] = 1;
	agent.Q[2][3][policyData[11]] = 1;
	agent.Q[2][4][policyData[12]] = 1;
	agent.Q[2][5][policyData[13]] = 1;
	agent.Q[2][6][policyData[14]] = 1;

	printTablesSideBySide(agent.Q, agent.N, "Q table of policy-" + to_string(policyNo), "N table of policy-" + to_string(policyNo));
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
	double startP = getCPUTime();
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

	sol.timeP = getCPUTime() - startP; // save the preprocessing time
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
	file << fixed << setprecision(2) << sol.policyNo << "\t" << filein << "\t" << sol.opt << "\t" << sol.timeT << "\t" << sol.iLB << "\t" << sol.LB
		<< "\t" << sol.Ncalls << "\t" << sol.Ncuts << "\t" << sol.NfailSl
		<< "\t" << sol.tLBinc << "\t" << sol.tInitSol
		<< "\t" << sol.tSlave << "\t" << sol.tRealSl
		<< "\t" << sol.tmis << "\t" << sol.tmisLR << "\t" << sol.tmisO << "\t" << sol.tmisR
		<< "\t" << sol.NmisLR << "\t" << sol.NmisLR_o << "\t" << sol.NmisO << "\t" << sol.NmisR << "\t" << sol.NmisFull
		<< "\t" << sol.DmisLR << "\t" << sol.DmisLR_o << "\t" << sol.DmisO << "\t" << sol.DmisR
		<< "\t" << sol.tCutCheck << "\t" << sol.NelCut << endl;
	file.close();
}


void writeVectorToFile(Agent& agent, const string& filenameQ, const string& filenameN) {
	ofstream outFile(filenameQ);

	if (outFile.is_open()) {
		size_t nbLevels = agent.Q.size();
		size_t nbState = agent.Q[1].size();
		size_t nbAct = agent.Q[0][0].size();

		// Write dimensions
		outFile << nbLevels << " " << nbState << " " << nbAct << endl;

		// Write data
		for (const auto& level : agent.Q) {
			for (const auto& state : level) {
				for (double act : state) {
					outFile << act << " ";
				}
				outFile << endl;
			}
			outFile << endl;
		}
	}

	outFile.close();


	ofstream outFileN(filenameN);

	if (outFileN.is_open()) {
		size_t nbLevels = agent.N.size();
		size_t nbState = agent.N[1].size();
		size_t nbAct = agent.N[0][0].size();

		// Write dimensions
		outFileN << nbLevels << " " << nbState << " " << nbAct << endl;

		// Write data
		for (const auto& level : agent.N) {
			for (const auto& state : level) {
				for (double act : state) {
					outFileN << act << " ";
				}
				outFileN << endl;
			}
			outFileN << endl;
		}
	}

	outFileN.close();

}


void readVectorFromFile(Agent& agent, const string& filenameQ, const string& filenameN) {
	int nbLevels, nbState, nbAct;

	ifstream infile(filenameQ);

	if (!infile.is_open()) {
		// create the empty ones
		int nbLevels = 3, nbState = 7, nbAct = 3;
		// Initialize the 3D vector with dimensions 5x7x2 filled with 0.0
		vector< vector< vector<double>>> Q(nbLevels, vector< vector<double>>(nbState, vector<double>(nbAct, 0.0)));
		// Modify Q[0] to have only 1 sub-vector of dimension 3
		Q[0] = vector< vector<double>>(1, vector<double>(nbAct, 10000.0));
		Q[1] = vector< vector<double>>(nbState, vector<double>(nbAct, 9500.0));
		Q[2] = vector< vector<double>>(nbState, vector<double>(2, 6000.0));
		//Q[3] = vector< vector<double>>(nbState, vector<double>(nbAct, 3500.0));

		// Initialize the 3D vector with dimensions 5x7x2 filled with 0
		vector< vector< vector<int>>> N(nbLevels, vector< vector<int>>(nbState, vector<int>(nbAct, 0)));
		// Modify Q[0] to have only 1 sub-vector of dimension 3
		N[0] = vector< vector<int>>(1, vector<int>(nbAct, 0));
		N[2] = vector< vector<int>>(nbState, vector<int>(2, 0));
		agent.Q = Q;
		agent.N = N;
		return;
	}

	infile >> nbLevels >> nbState >> nbAct;
	vector< vector< vector<double>>> table(nbLevels);
	for (int i = 0; i < nbLevels; ++i) {
		int stateSize = (i == 0) ? 1 : nbState;
		int actSize = (i == nbLevels - 1) ? 2 : nbAct;
		table[i] = vector< vector<double>>(stateSize, vector<double>(actSize));

		for (int j = 0; j < stateSize; ++j) {
			for (int k = 0; k < actSize; ++k) {
				infile >> table[i][j][k];
			}
		}
	}
	agent.Q = table;
	infile.close();


	ifstream infileN(filenameN);
	infileN >> nbLevels >> nbState >> nbAct;
	vector< vector< vector<int>>> N(nbLevels);
	for (int i = 0; i < nbLevels; ++i) {
		int stateSize = (i == 0) ? 1 : nbState;
		int actSize = (i == nbLevels - 1) ? 2 : nbAct;
		N[i] = vector< vector<int>>(stateSize, vector<int>(actSize));

		for (int j = 0; j < stateSize; ++j) {
			for (int k = 0; k < actSize; ++k) {
				infileN >> N[i][j][k];
			}
		}
	}
	agent.N = N;
	infileN.close();

}

void printTablesSideBySide(const vector<vector<vector<double>>>& Q, const vector<vector<vector<int>>>& N, const string& name1, const string& name2) {

	for (size_t i = 0; i < Q.size(); ++i) {
		cout << "Q[" << i << "]:" << endl;
		for (size_t j = 0; j < Q[i].size(); ++j) {
			cout << "  Q[" << i << "][" << j << "]: ";
			for (size_t k = 0; k < Q[i][j].size(); ++k) {
				cout << fixed << setprecision(2) << setw(7) << Q[i][j][k] << " ";
			}
			if (i == Q.size() - 1) cout << setw(8) << "   ";
			cout << "  N[" << i << "][" << j << "]: ";
			for (size_t k = 0; k < N[i][j].size(); ++k) {
				cout << setw(2) << N[i][j][k] << " ";
			}
			cout << endl;
		}
	}
}

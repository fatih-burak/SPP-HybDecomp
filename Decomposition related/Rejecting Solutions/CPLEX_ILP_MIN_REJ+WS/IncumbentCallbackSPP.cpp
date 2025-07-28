using namespace std;
#include "IncumbentCallbackSPP.h"
#include "model.h"
#include "helper_functions.h"
#include "slave_functions.h"


IncumbentCallbackSPP::IncumbentCallbackSPP(IloEnv env, Instance& instS, Solution& solS, const vector<vector<IloIntVar>>& xS, const IloIntVar zS, const vector<vector<bool>>& NPsS)
    : IloCplex::IncumbentCallbackI(env), inst(instS), sol(solS), x(xS), z(zS), NPs(NPsS)
{
    // Constructor body (if needed)
}

void IncumbentCallbackSPP::main() {
    try {
		int optz = sol.iLB;
		//cout << getCPUTime() - sol.startTime << endl;
		if (getCPUTime() - sol.startTime + 1 > inst.timeLimit || sol.opt) {
			if (sol.opt == 0) cout << "Time limit reached: abort... " << endl;
			abort(); return;
		}

		// Get the best current lower bound from CPLEX
		double currentBestLowerBound = getBestObjValue();
		cout << "Best LB: " << currentBestLowerBound;
		//get current height
		int H = getValue(z);
		cout << " vs trying UB = " << H << endl;

		if (sol.tInitSol == -1) sol.tInitSol = getCPUTime() - sol.startTime;
		sol.Ncalls += 1; // increase the number of calls by 1
		cout << "Ncalls " << sol.Ncalls << " calls " << endl;

        int n2 = inst.items2.size(); //nb of items
		vector<int> getStart(n2);
		for (int j = 0; j < n2; j++) {
			for (int i = 0; i <= inst.W - inst.items2[j][0]; i++) {
				if (NPs[j][i] && ceil(getValue(x[i][j]) - EPSILON) > 0) getStart[j] = i;
			}
		}

		vector<vector<int> > items_slave(n2);
		for (int j = 0; j < n2; j++) items_slave[j] = { inst.items2[j][0],inst.items2[j][1], getStart[j] }; //width, height, x-coord

		vector<vector<int>> bins(inst.W); //create a vector which stores the vertical item pieces packed in each column (bin)
		for (int j = 0; j < n2; j++) {
			for (int k = getStart[j]; k < getStart[j] + items_slave[j][0]; k++) bins[k].push_back(j);
		}
		//print_master_sol(bins); // print the master solution per bin

		// CALL THE SLAVE
		vector<vector<int>> current_sol;
		double startSlave = getCPUTime(); //to measure the time spent in the slave
		int status = slave(items_slave, bins, H, current_sol, 0);
		cout << fixed << setprecision(2) << "H = "<< H <<" and slave status = " << status << " and time in slave = " << getCPUTime() - startSlave << endl;
		if (status == 1) { // feasible 
			if (H > optz) {
				sol.higherCuts += 1;
				sol.higherFeas += 1;
				sol.t_higherFeas += getCPUTime() - startSlave;
				sol.t_higherCuts += getCPUTime() - startSlave;
			}

			cout << "Slave feasible." << endl; //but since minimization, you cannot end the search.
			//update UB
			if (sol.UB > H) sol.UB = H;
			if (sol.UB == ceil(currentBestLowerBound-EPSILON) || sol.UB == sol.iLB) {
				sol.opt = 1;
				// Stop the CPLEX optimization process
				abort(); return;
			}
		}
		else if (status == 0) { //time limit reached
			if (H > optz) {
				sol.higherCuts += 1;
				sol.higherSkip += 1;
				sol.t_higherSkip += getCPUTime() - startSlave;
				sol.t_higherCuts += getCPUTime() - startSlave;
			}

			cout << "Time limit in slave: skipping the solution..." << endl;
			//sol.cuts.insert({ one_cut });
			sol.NfailSl++; //increase the nb fail in slave (This number is not included in Ncuts)
			cout << "NfailSl added " << sol.NfailSl << " cuts..." << endl;
			cout << "-------------------------------------------------------------------" << endl;
			reject();
		}
		else if (status == 3) { //infeasible
			if (H > optz) {
				sol.higherCuts += 1;
				sol.higherInfeas += 1;
				sol.t_higherInfeas += getCPUTime() - startSlave;
				sol.t_higherCuts += getCPUTime() - startSlave;
			}
			cout << "Slave infeasible..." << endl;
			sol.Ncuts++;
			reject();
		}
		cout << "Ncalls = " << sol.Ncalls << ": infeasible vs failed = " << sol.Ncuts << ", " << sol.NfailSl << endl;
		cout << "-------------------------------------------------------------------" << endl;
		// END SLAVE //
		sol.tSlave += getCPUTime() - startSlave; //to measure the time spent in the slave
		cout << "Time Slave=" << getCPUTime() - startSlave << endl;
		cout << "-------------------------------------------------------------------" << endl;
		// Terminate Gurobi when time limit in callback is reached.
		cout << "Time info = " << getCPUTime() - sol.startTime << endl;
		if (getCPUTime() - sol.startTime + 1 > inst.timeLimit) {
			cout << "Time limit reached: abort..." << endl;
			abort(); return;
		}
    }
    catch (IloException& e) {
        cerr << "Error in callback: " << e << endl;
    }
}
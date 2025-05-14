#include "RL.h"
#include <random>

void train(Instance& inst, Solution& sol, Agent& agent) {
	double startT = getCPUTime();
	// Q parameters
	//double gamma = 1; //discount factor
	//double alpha = 0.1; // learning rate
	double epsilon = 0.1; // greedy move prob

	cout << "Agent is trying H = " << sol.LB + sol.prepacked_height << " actually trying H (w/o preprocessed height)= " << sol.LB << endl; //show the tried height together with the prepacked height
	int idx = 0; //which level you are at (kind of which state you are at)
	int idx2 = 0; //which model-state you are at
	int a; //action
	//double r; //reward

	while (idx < agent.Q.size() && sol.opt == 0) {
		cout << "\n==========================================================================================================" << endl;
		cout << "\n----------------------------------------------------------------------------------------------------------" << endl;
		cout << "Agent is now at:" << idx << " " << idx2 << endl;
		//FIND THE ACTION
		random_device rd;  // Obtain a random number from hardware
		mt19937 gen(rd()); // Seed the generator
		uniform_real_distribution<> dis(0.0, 1.0);
		double random_number = dis(gen);  // Generate a uniformly distributed random number between 0 and 1
		//cout << "Random number = " << random_number << endl;
		if (0 && random_number < epsilon) { //deterministic behaviour
			//EXPLORE: take a random action
			cout << "Agent takes a random action..." << endl;
			a = rand() % agent.Q[idx][idx2].size(); // 0: go to ILP, 1: go to CP, 2: go to both
		}
		else {
			// Find the maximum element in agent.Q[idx][idx2]
			auto max_it = max_element(agent.Q[idx][idx2].begin(), agent.Q[idx][idx2].end());
			// Calculate the index of the maximum element
			a = distance(agent.Q[idx][idx2].begin(), max_it);
		}
		cout << "Selected action = " << a << endl;
		cout << "----------------------------------------------------------------------------------------------------------" << endl;
		int idxb = idx, idx2b = idx2; //SAVE THE STATE YOU ARE IN
		int prevNbSol = sol.Ncalls;
		//UPDATE N, PERFORM THE ACTION AND UPDATE THE STATE
		agent.N[idx][idx2][a] += 1;
		double startRL = getCPUTime();
		RL(inst, sol, idx, idx2, a);
		double actualTime = getCPUTime() - startRL;
		sol.totalTriedTimes += actualTime;
		//OBSERVE THE REWARD
		if (sol.opt == 0 && sol.isLBincreased) {
			cout << "----------------------------------------------------------------------------------------------------------" << endl;
			cout << "H is increased!" << endl;
			cout << "----------------------------------------------------------------------------------------------------------" << endl;
			cout << "Time it takes to solve the model = " << round((actualTime) * 100.0) / 100.0 << endl;
			double time2Solve = getCPUTime() - startT;
			//cout << "Observed (total) time for the reward = " << round(time2Solve * 100.0) / 100.0 << endl;
			//r = 9000 * exp(-0.01 * time2Solve) + 1000;
			//cout << "Observed reward = " << round(r * 100.0) / 100.0 << endl;
		}
		else if (sol.opt == 1) {
			cout << "Time it takes to solve the model = " << round(actualTime * 100.0) / 100.0 << endl;
			double time2Solve = getCPUTime() - startT;
			//cout << "Observed (total) time for the reward = " << round(time2Solve * 100.0) / 100.0 << endl;
			//r = 9000 * exp(-0.01 * time2Solve) + 1000;
			//cout << "Observed reward = " << round(r * 100.0) / 100.0 << endl;
		}
		else if (sol.opt == 0) {
			// UPTADE THE agent.Q-VAL
			cout << "Agent could not solve the problem..." << endl;
			//cout << "Agent could not solve the problem and will be rewarded for the solutions it found..." << endl;
			int nbNewSolsFound = sol.Ncalls - prevNbSol;
			//r = nbNewSolsFound;
			//cout << "NB sols found =  " << r << endl;
		}

		/*
		//UPDATE THE agent.Q TABLE
		cout << "Updating: Q[" << idxb << "][" << idx2b << "][" << a << "]" << endl;
		if (sol.opt || idx == agent.Q.size() || (sol.opt == 0 && sol.isLBincreased)) { //if opt is found or agent is in the last state or (opt could not found but LB is increase)
			cout << "Old Q[" << idxb << "][" << idx2b << "][" << a << "] val = " << agent.Q[idxb][idx2b][a] << endl;
			agent.Q[idxb][idx2b][a] = agent.Q[idxb][idx2b][a] * (1 - alpha) + alpha * r;
			cout << "New Q[" << idxb << "][" << idx2b << "][" << a << "] val = " << agent.Q[idxb][idx2b][a] << endl;
		}
		else {
			double maximum_element = *max_element(agent.Q[idx][idx2].begin(), agent.Q[idx][idx2].end());
			cout << "maximum q-val in the row:" << idx << " " << idx2 << " is " << maximum_element << endl;
			cout << "Old Q[" << idxb << "][" << idx2b << "][" << a << "] val = " << agent.Q[idxb][idx2b][a] << endl;
			agent.Q[idxb][idx2b][a] = agent.Q[idxb][idx2b][a] * (1 - alpha) + alpha * (r + gamma * maximum_element);
			cout << "New Q[" << idxb << "][" << idx2b << "][" << a << "] val = " << agent.Q[idxb][idx2b][a] << endl;
		}

		// Print the values of Q nicely
		printTablesSideBySide(agent.Q, agent.N, "Q", "N");
		*/
		cout << "Is it solved? = " << sol.opt << endl;

		if (sol.opt == 0 && sol.isLBincreased) {
			cout << "----------------------------------------------------------------------------------------------------------" << endl;
			cout << "Agent is going back to the START state." << endl;
			cout << "----------------------------------------------------------------------------------------------------------" << endl;
			cout << "Agent is trying H = " << sol.LB + sol.prepacked_height << " actually trying H (w/o preprocessed height)= " << sol.LB << endl; //show the tried height together with the prepacked height
			// return back to state 0 with increased H
			//reset values
			idx = 0; idx2 = 0;
			startT = getCPUTime(); //START THE CLOCK
			sol.isLBincreased = false;
		}
	}
}

void RL(Instance& inst, Solution& sol, int& idx, int& idx2, int& a) {
	//MOVE TO THE NEXT LEVEL
	idx++;
	int prevNbSol = sol.Ncalls;
	
	// SELECT THE TIME LIMIT
	if (idx == 1) {
		inst.timeLimit = 6; // if at level 1
	}
	else if (idx == 2) {
		inst.timeLimit = 60;  // if at level 2
	}
	else if (idx == 3) {
		//inst.timeLimit = max(inst.totalTimeRem - sol.totalTriedTimes, 0.001);  //if it is the last level (level 4), give the full amount of remaining time
		inst.timeLimit = 1200-66;  //if it is the last level (level 4), give the full amount of remaining time
	}

	createStory(idx, idx2, a, inst.timeLimit); //this is based on where the agent was before taking the action.

	if (sol.totalTriedTimes > inst.totalTimeRem) { cout << "Time limit reached!" <<endl; return; } //time limit reached

	// SELECT THE METHOD
	int nbSolBM = 0, nbSolCP = 0;
	if (a == 0) {
		cout << "agent is trying ILP for " << inst.timeLimit << " seconds." << endl;
		MP_ILP(inst, sol); // find the BM solution
	}
	else if (a == 1) {
		cout << "agent is trying CP for " << inst.timeLimit << " seconds." << endl;
		MP_CP(inst, sol); // find the CPM solution 
	}
	else if (a == 2) {
		inst.timeLimit /= 2;
		double random_number = (double)rand() / RAND_MAX; // get the random number
		if (random_number < 0.5) {
			cout << "agent is trying ILP first for " << inst.timeLimit << " seconds." << endl;
			MP_ILP(inst, sol); // find the BM solution
			nbSolBM = sol.Ncalls - prevNbSol;
			cout << "\tfound " << nbSolBM << " solutions." << endl;
			if (sol.opt || sol.isLBincreased) return; // return either when optimal sol is found or LB is increased
			cout << "agent is now trying CP for " << inst.timeLimit << " seconds." << endl;
			MP_CP(inst, sol); // find the CPM solution
			nbSolCP = sol.Ncalls - prevNbSol - nbSolBM;
			cout << "\tfound " << nbSolCP << " solutions." << endl;
		}
		else {
			cout << "agent is trying CP first for " << inst.timeLimit << " seconds." << endl;
			MP_CP(inst, sol); // find the CPM solution
			nbSolCP = sol.Ncalls - prevNbSol;
			cout << "\tfound " << nbSolCP << " solutions." << endl;
			if (sol.opt || sol.isLBincreased) return; // return either when optimal sol is found or LB is increased
			cout << "agent is now trying ILP for " << inst.timeLimit << " seconds." << endl;
			MP_ILP(inst, sol); // find the BM solution
			nbSolBM = sol.Ncalls - prevNbSol - nbSolCP;
			cout << "\tfound " << nbSolBM << " solutions." << endl;
		}

	}
	cout << "\tfound " << sol.Ncalls - prevNbSol << " solutions (in total)." << endl;
	if (sol.isLBincreased) return; // return either when optimal sol is found or LB is increased

	//FIND NEXT STATE (where you ended up after solving the model(s)): Update idx2
	if (a == 2) { //if you are at the combined state
		if (nbSolBM - nbSolCP > 0) idx2 = a * 2;   //(-1,-,-)
		else if (nbSolBM - nbSolCP == 0) idx2 = a * 2 + 1; // (-,0,-)
		else idx2 = a * 2 + 2; // (-,-,1)
	}
	else {
		if (sol.Ncalls - prevNbSol > 0) idx2 = a * 2 + 1;
		else idx2 = a * 2;
	}
}

void createStory(int& idx, int& idx2, int& a, double& time_Val) {
	cout << "----------------------------------------------------------------------------------------------------------" << endl;
	cout << "Agent is at ";
	if (idx - 1 == 0) cout << "the START ";
	else {
		cout << "level-" << idx - 1 << " state-" << idx2 << ": found ";
		if (idx2 == 0 || idx2==2) cout << "0 solutions and ";
		if (idx2 == 1 || idx2 == 3) cout << "more than or equal to 1 solutions and ";
	}
	cout << "chooses action ";
	if (a == 0) cout << "ILP for " << time_Val << " seconds." << endl;
	else if (a == 1) cout << "CP for " << time_Val << " seconds." << endl;
	else if (a == 2) cout << "both ILP then CP for " << time_Val / 2 << " seconds." << endl;


	cout << "----------------------------------------------------------------------------------------------------------" << endl;
}
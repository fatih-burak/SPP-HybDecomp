#include "main.h"
#include <thread>

int main(int argc, char **argv){          
    	// Read input and output paths
	int policyNo = stoi(argv[1]);	
	string pathPolicyFile = argv[2];
	string path = argv[3];	
	string filein = argv[4];
	string pathAndFileout = argv[5];

	Solution sol(filein);
	sol.policyNo = policyNo ;
	cout << filein << endl;	sol.name = filein;

	Agent agent;
	readPolicyVectorFromFile(agent, policyNo, pathPolicyFile);

    	// initialize the input variables from a file
	double startP = getCPUTime();
    	Instance inst = readInstance(path + filein, sol);
	sol.timeP = getCPUTime() - startP;
	inst.print();
	cout << "Size: " << sol.saved_cuts.size() << endl; // Should print 0

	inst.totalTimeRem= 1205;
	double start = getCPUTime();
	sol.startT = start;

    	// find the LB solution
	LB(inst, sol);
	cout << "Initial LB = " << sol.iLB << endl;

	// find the model solution
	train(inst, sol, agent);

	//add the prepacked height
	sol.LB += sol.prepacked_height;
	sol.iLB += sol.prepacked_height;
	sol.timeT = getCPUTime() - start;
	printInfo(pathAndFileout, sol, filein);
	printTablesSideBySide(agent.Q, agent.N, "Q", "N");
}



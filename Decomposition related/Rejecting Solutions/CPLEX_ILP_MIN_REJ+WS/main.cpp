#include "main.h"
#include <thread>

int main(int argc, char **argv){          
    // Read input and output paths
	string path = argv[1];	
	string filein = argv[2];
	string pathAndFileout = argv[3];
	Solution sol(filein);

	double start = getCPUTime();

	cout << filein << endl;	sol.name = filein;
    // initialize the input variables from a file
	double startP = getCPUTime();
    Instance inst = readInstance(path + filein, sol);
	inst.timeLimit = 1200; 
	sol.timeP = getCPUTime() - startP;
	inst.print();

    // find the LB solution
	LB(inst, sol);
	sol.LB = sol.iLB;
	BEST_FIT(inst, sol);
	cout << "Initial LB = " << sol.iLB << " and prepacked = "<< sol.prepacked_height<< endl;

	// find the model solution
	ILP_CPLEX_MIN_REJ(inst, sol);
	//add the prepacked height
	sol.LB += sol.prepacked_height;
	sol.iLB += sol.prepacked_height;
	sol.timeT = getCPUTime() - start;
	printInfo(pathAndFileout, sol, filein);
}



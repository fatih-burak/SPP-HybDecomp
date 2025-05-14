#include "main.h"
#include <thread>

int main(int argc, char **argv){          
    // Read input and output paths
	string path = argv[1];	
	string filein = argv[2];
	string pathAndFileout = argv[3];
	Solution sol;

	cout << filein << endl;	sol.name = filein;
    // initialize the input variables from a file
	double startP = getCPUTime();
    	Instance inst = readInstance(path + filein, sol);
	inst.timeLimit = 1200; 
	sol.timeP = getCPUTime() - startP;
	inst.print();

    // find the LB solution
	LB(inst, sol);
	cout << "Initial LB = " << sol.iLB << endl;
	double start = getCPUTime();
	// find the model solution
	MP_ILP_DEL_LCBC(inst, sol);
	//add the prepacked height
	sol.LB += sol.prepacked_height;
	sol.iLB += sol.prepacked_height;
	sol.timeT = getCPUTime() - start;
	printInfo(pathAndFileout, sol, filein);
}



#include "main.h"
#include <thread>

int main(int argc, char **argv){          
	// Read input and output paths
	string path = argv[1];	
	string filein = argv[2];
	string pathAndFileout = argv[3];
	Solution sol; 
	
    // initialize the input variables from a file
	double startP = getCPUTime();
    	Instance inst = readInstance(path + filein, sol);
	inst.name = filein;
	inst.timeLimit = 1200*2; 
	sol.timeP = getCPUTime() - startP;
	//inst.print();

     // find the LB solution
	LB(inst, sol);
	sol.LB = sol.iLB;
	cout << "Initial LB = " << sol.iLB << " and prepacked = "<< sol.prepacked_height<< endl;

     // find the solution
	double start = getCPUTime();
	CPM(inst, sol);

     // add the prepacked height
	sol.LB += sol.prepacked_height;
	sol.iLB += sol.prepacked_height;

	sol.timeT = getCPUTime() - start;
	printInfo(pathAndFileout, sol, filein);
	//printSInfo(sol, inst);
}



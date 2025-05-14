#include "main.h"
#include <thread>

int main(int argc, char **argv){          
    // Read input and output paths
	string path = argv[2];	
	string filein = argv[3];
	string pathAndFileout = argv[4];
	string path_UB = argv[1];
	Solution sol; 

	double start = getCPUTime();	

    // initialize the input variables from a file
	double startP = getCPUTime();
    	Instance inst = readInstance(path + filein, sol);
	inst.timeLimit = 1200; 
	sol.timeP = getCPUTime() - startP;
	inst.print();
    // find the LB solution
	LB(inst, sol);
	sol.LB = sol.iLB;
    // read the UB
	filein.erase(filein.size() - 4);
	cout << filein << endl;
	sol.UB = getValueForInstance(path_UB, filein) - sol.prepacked_height; 
	cout << "Initial LB = " << sol.iLB << " and UB = " << sol.UB << " and prepacked = "<< sol.prepacked_height<< endl;
    // find the solution
	CC_ILP(inst, sol);
	
    //add the prepacked height
	sol.LB += sol.prepacked_height;
	sol.iLB += sol.prepacked_height;
	sol.UB += sol.prepacked_height; 

	sol.timeT += getCPUTime() - start;
	printInfo(pathAndFileout, sol, filein);
}



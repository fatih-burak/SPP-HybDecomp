#include "main.h"
#include <thread>

int main(int argc, char **argv){          
    // Read input and output paths
	string path = argv[1];	
	string filein = argv[2];
	string pathAndFileout = argv[3];
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
    // get the UB
	filein.erase(filein.size() - 4);
	cout << filein << endl;
	BEST_FIT(inst, sol);
	cout << "Initial LB = " << sol.iLB << " and UB = " << sol.UB << " and prepacked = "<< sol.prepacked_height<< endl;
    // find the solution
	CC_CP(inst, sol);
	
    //add the prepacked height
	sol.LB += sol.prepacked_height;
	sol.iLB += sol.prepacked_height;
	sol.UB += sol.prepacked_height; 

	sol.timeT += getCPUTime() - start;
	printInfo(pathAndFileout, sol, filein);
}



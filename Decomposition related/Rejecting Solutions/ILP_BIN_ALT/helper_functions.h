#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

using namespace std;
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream> 
#include <time.h>
//#include <sys/time.h>
#include <ilcp/cp.h>
#include <set>
#include <math.h> 
#include <iostream>

const double EPSILON = 0.00001; // small constant
const int M = 10000;

struct Instance
{
	int n; 						// number of items
	int W, Wo; 					// number of bins
	vector<vector<int>> items;	// items
	vector<vector<int>> items2;	// binary itemsQ
	vector<vector<bool> > NPs; // normal patterns
	double timeLimit = 1200;
	void print();
};

struct Solution
{
	double startTime;
	//General
	string name = "instance";
	int opt = 0, prepacked_height = 0, Ncuts = 0, Ncalls = 0, NfailSl = 0;
	double timeT = 0.0, tSlave = 0.0, tRealSl = 0.0, tLBinc = 0.0, tFeasSl = 0.0, timeP = 0.0, tInitSol = -1.0;
	//slave preprocess
	int NSlInc_w = 0, NSDec_W = 0;

	//LB related
	int iLB = 0, LB = 0;
	int vLB1 = 0, vLB2 = 0, vLB3 = 0;
	double tLB = 0.0, tLB1 = 0.0, tLB2 = 0.0, tLB3 = 0.0;

	vector<vector<int>> assignedBin;
	vector<vector<int>> solution;
	set<set<vector<int>>> cuts;

	vector<vector<int>> all_sols; //!
	vector<double> percentages; //!
	double perc_sim = -1.0; //!

	Solution(string name_val) : name(name_val) {}

};

double getCPUTime();

Instance readInstance(const string& filename, Solution& sol);
int CP(const Instance& inst, const int& heightL, vector<vector<int> >& assignedBin);
void printInfo(const string& pathAndFileout, const Solution& sol, const string& filein);

//SORTING FUNCTIONS
bool sortItems(const vector<int>& v1, const vector<int>& v2);
bool sortItemsArea(const vector<int>& v1, const vector<int>& v2);
bool sortWidthSL(const vector<int>& v1, const vector<int>& v2);
bool sortItemsStartPosLS(const std::vector<int>& v1, const std::vector<int>& v2);
bool sortItemsStartPosSL(const std::vector<int>& v1, const std::vector<int>& v2);

#endif 

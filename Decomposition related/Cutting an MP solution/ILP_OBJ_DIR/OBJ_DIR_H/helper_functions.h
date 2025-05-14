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
#include <ilcp/cp.h>
#include "gurobi_c++.h"
#include <set>
#include <math.h> 
#include <iostream>

const double EPSILON = 0.00001; // small constant
const int M = 1000000;			// big constant

struct Instance
{
	int n; 						// number of items
	int W, Wo; 					// number of bins
	vector<vector<int>> items;	// items
	vector<vector<int>> items2;	// binary itemsQ
	vector<vector<bool> > NPs; // normal patterns
	double timeLimit = 1200;
	double totalTimeRem = 1200;
	void print();
};

struct Solution
{
	//General
	string name = "instance";
	int opt = 0, prepacked_height = 0, Ncuts = 0, Ncalls = 0, NfailSl = 0;
	double timeT = 0.0, tSlave = 0.0, tRealSl = 0.0, tLBinc = 0.0, timeP = 0.0, tInitSol = -1.0;
	//slave preprocess
	int NSlInc_w = 0, NSDec_W = 0;

	//LB related
	int iLB = 0, LB = 0;
	int vLB1 = 0, vLB2 = 0, vLB3 = 0;
	double tLB = 0.0, tLB1 = 0.0, tLB2 = 0.0, tLB3 = 0.0;

	//MIS related
	int NmisLR = 0, NmisLR_o = 0, NmisO = 0, NmisR = 0, NmisFull = 0; //how many times succes. cuts are added
	int DmisLR = 0, DmisLR_o = 0, DmisO = 0, DmisR = 0; //how many deleted
	double tmis = 0.0, tmisLR = 0.0, tmisO = 0.0, tmisR = 0.0; //time spent

	//eliminated cuts due to copy or super set
	int NelCut = 0; //nb of times found
	double tCutCheck = 0.0; //time spent checking

	vector<vector<int>> assignedBin;
	vector<vector<int>> solution;
	set<set<vector<int>>> cuts;
	vector<vector<int>> all_sols; //!
	vector<double> percentages; //!
	double perc_sim = -1.0; //!
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

#ifndef SLAVE_FUNCTIONS_H 
#define SLAVE_FUNCTIONS_H
#include <random>
#include "helper_functions.h"
#include "model.h"
//#include "gurobi_c++.h"

int slave(const vector<vector<int>>& items_slave, const vector<vector<int>>& bins, const int& H, vector<vector<int>>& current_sol, int type); //type ensures different behaviour for slave. 0: main slave, 1: mis, -1: preprocess slave

int slave_CP(const vector<vector<int>>& items_slave, const vector<vector<int>>& bins, const int& H, vector<vector<int>>& current_sol, int type, string instanceName); //type ensures different behaviour for slave. 0: main slave, 1: mis, -1: preprocess slave
int slave_ILP(const vector<vector<int>>& items_slave_org, const vector<vector<int>>& bins_org, const int& H, vector<vector<int>>& current_sol, int type, string instanceName);

void preprocess_slave(vector<vector<int>>& items_slave, vector<vector<int>>& bins);
void preprocess_slave2(vector<vector<int>>& items_slave, vector<vector<int>>& bins);

void print_master_sol(const vector<vector<int>>& bins); //prints the solution of the master per bin
vector<bool> MIS_LR(const vector<vector<int>>& items_slave, const vector<vector<int>>& bins, const int& H); //MIS: left, right delete
vector<bool> MIS_O(const vector<vector<int>>& items_slave, const vector<vector<int>>& bins, const int& H, vector<bool>& deleted); //MIS: ordered
vector<bool> MIS_R(const vector<vector<int>>& items_slave, const vector<vector<int>>& bins, const int& H, vector<bool>& deleted, int seed_value); //MIS: random
vector<vector<int>> LIFTING(Instance& inst, const vector<vector<int>>& items_slave, vector<bool>& deleted);
bool cut_check(const vector<vector<int>>& items_slave, vector<bool>& deleted, Solution& sol);
void re_add_cuts(const vector<vector<int>>& items2, const int& W, const int& H, Solution& sol);

#endif

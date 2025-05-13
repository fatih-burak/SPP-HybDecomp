#ifndef SLAVE_FUNCTIONS_H 
#define SLAVE_FUNCTIONS_H

#include "helper_functions.h"
#include "model.h"
#include <random>

int slave(const vector<vector<int>>& items_slave, const vector<vector<int>>& bins, const int& H, vector<vector<int>>& current_sol, int type);
void preprocess_slave(vector<vector<int>>& items_slave, vector<vector<int>>& bins);

void print_master_sol(const vector<vector<int>>& bins); //prints the solution of the master per bin
vector<bool> MIS_LR(const vector<vector<int>>& items_slave, const vector<vector<int>>& bins, const int& H); //MIS: left, right delete
vector<bool> MIS_O(const vector<vector<int>>& items_slave, const vector<vector<int>>& bins, const int& H, vector<bool>& deleted); //MIS: ordered
vector<bool> MIS_R(const vector<vector<int>>& items_slave, const vector<vector<int>>& bins, const int& H, vector<bool>& deleted, const int& seed_value); //MIS: random
vector<vector<int>> LIFTING(Instance& inst, const vector<vector<int>>& items_slave, vector<bool>& deleted);
bool cut_check(const vector<vector<int>>& items_slave, vector<bool>& deleted, Solution& sol);
void re_add_cuts(const vector<vector<int>>& items2, const int& W, const int& H, Solution& sol);

#endif 
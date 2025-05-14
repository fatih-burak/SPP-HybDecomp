#ifndef SLAVE_FUNCTIONS_H 
#define SLAVE_FUNCTIONS_H

#include "helper_functions.h"
#include "model.h"
#include <random>

int slave(const vector<vector<int>>& items_slave, const vector<vector<int>>& bins, const int& H, vector<vector<int>>& current_sol, int type); //type ensures different behaviour for slave. 0: main slave, 1: mis, -1: preprocess slave
void preprocess_slave2(vector<vector<int>>& items_slave, vector<vector<int>>& bins);
void print_master_sol(const vector<vector<int>>& bins); //prints the solution of the master per bin

bool cut_check(const vector<vector<int>>& items_slave, vector<bool>& deleted, Solution& sol);
void re_add_cuts(const vector<vector<int>>& items2, const int& W, const int& H, Solution& sol);

#endif 
#ifndef SLAVE_FUNCTIONS_H 
#define SLAVE_FUNCTIONS_H

#include "helper_functions.h"
#include "model.h"
#include <random>

int slave(const vector<vector<int>>& items_slave, const vector<vector<int>>& bins, const int& H, vector<vector<int>>& current_sol, int type);
void preprocess_slave(vector<vector<int>>& items_slave, vector<vector<int>>& bins);

#endif 
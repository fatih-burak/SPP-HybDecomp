#ifndef SLAVE_FUNCTIONS_H 
#define SLAVE_FUNCTIONS_H
#include <random>
#include "helper_functions.h"
#include "model.h"

int slave(const vector<vector<int>>& items_slave, const vector<vector<int>>& bins, const int& H, vector<vector<int>>& current_sol, int type); //type ensures different behaviour for slave. 0: main slave, 1: mis, -1: preprocess slave

#endif

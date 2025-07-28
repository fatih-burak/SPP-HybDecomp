#ifndef MODEL_H
#define MODEL_H

using namespace std;
#include <iostream>


#include "gurobi_c++.h"
#include "helper_functions.h"
#include "slave_functions.h"

void DC_ILP(Instance& inst, Solution& sol);
void BEST_FIT(Instance& inst, Solution& sol);

#endif 


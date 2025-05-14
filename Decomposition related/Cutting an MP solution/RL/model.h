#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include "gurobi_c++.h"
#include "helper_functions.h"
#include "slave_functions.h"
using namespace std;

void MP_ILP(Instance& inst, Solution& sol);
void MP_CP(Instance& inst, Solution& sol);

#endif 

